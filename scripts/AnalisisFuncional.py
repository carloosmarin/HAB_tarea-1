#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Tarea 1: Análisis Funcional (COX4I2, ND1, ATP6)

Este script realiza un análisis funcional de una lista pequeña de genes (por defecto,
lee COX4I2, ND1 y ATP6 desde data/genes_input.txt). Está pensado para listas cortas,
donde el enfoque recomendado es el análisis de sobre-representación (ORA).

Métodos y bases de datos utilizados (solo los vistos en clase/conversación):
- MyGene.py (paquete: mygene) para normalizar símbolos humanos y convertir a UniProt/Ensembl/Entrez.
  Fuente agregada de anotaciones (MyGene.info).
- Biopython/Entrez (opcional) para recuperar resúmenes/metadata desde NCBI Gene (requiere email).
- GOATOOLS para ORA de Gene Ontology (GO), usando:
    * Ontología: go-basic.obo (descarga desde Gene Ontology)
    * Asociaciones: goa_human.gaf (UniProtKB→GO; descarga oficial)
  Estadística: Test exacto de Fisher + corrección múltiple BH (FDR).
- Enrichr (API) para ORA contra bibliotecas amplias (GO Biological Process, Reactome, KEGG_2019_Human).
  Devuelve p-value, Z-score y combined score (útil para priorizar).
- STRINGdb (API) para ORA de GO (y otras categorías) con respaldo de red de PPIs de STRING.

Flujo:
1) Leer genes de entrada (símbolos).
2) Normalizar/validar (humano) y mapear a UniProt/Ensembl/Entrez con mygene.
3) (Opcional) Descargar resúmenes de NCBI Gene con Biopython/Entrez.
4) Descargar go-basic.obo y goa_human.gaf; construir universo (background) y asociaciones (UniProt→GO).
5) ORA de GO con GOATOOLS (propagate_counts=True; FDR BH).
6) ORA con Enrichr (API) en bibliotecas seleccionadas (por defecto: GO BP y Reactome).
7) ORA con STRING (API) y guardar resultados.
8) Guardar todos los outputs en 'results/'.

CLI:
    python scripts/tu_script.py --input data/genes_input.txt --outdir results \
        --email_topt TU_EMAIL@DOMINIO --string_id TU_EMAIL@DOMINIO

Requisitos (requirements.txt sugerido):
    mygene
    biopython
    goatools
    requests
    pandas

Nota: GSEA (gseapy) NO se usa aquí porque con 3 genes no es adecuado un análisis por ranking.
"""

import argparse
import os
import sys
import time
import gzip
import shutil
import logging
import json
from typing import List, Dict, Any, Tuple

import requests
import pandas as pd
import mygene

# Biopython (Entrez) es opcional: solo si se pasa --email_topt
from Bio import Entrez

# GOATOOLS
from goatools import obo_parser
from goatools.associations import read_gaf
from goatools.go_enrichment import GOEnrichmentStudy


# ------------------------ Utilidades generales ------------------------ #

def setup_logging():
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s [%(levelname)s] %(message)s",
        datefmt="%H:%M:%S"
    )


def ensure_dir(path: str):
    os.makedirs(path, exist_ok=True)


def read_gene_list(input_path: str) -> List[str]:
    genes = []
    with open(input_path, "r", encoding="utf-8") as f:
        for line in f:
            g = line.strip()
            if g and not g.startswith("#"):
                genes.append(g)
    return genes


def save_df(df: pd.DataFrame, path: str):
    if df is None or df.empty:
        logging.warning(f"No hay datos para guardar en {path}")
        return
    df.to_csv(path, sep="\t", index=False)
    logging.info(f"Guardado: {path}")


# ------------------------ Paso 1: MyGene (normalización / mapping) ------------------------ #

def mygene_mapping(genes: List[str], species: str = "human") -> pd.DataFrame:
    """
    Normaliza símbolos humanos y obtiene múltiples IDs (UniProt, Ensembl, Entrez).
    - Usa scopes amplios ('symbol,alias') y species='human' para resolver alias (ND1 -> MT-ND1).
    - Unifica posibles múltiples UniProt (Swiss-Prot y/o TrEMBL).
    """
    mg = mygene.MyGeneInfo()
    logging.info(f"Consultando MyGene para {len(genes)} genes...")
    res = mg.querymany(
        genes,
        scopes="symbol,alias",
        fields="symbol,name,entrezgene,ensembl.gene,uniprot",
        species=species,
        as_dataframe=False
    )

    rows = []
    for r in res:
        q = r.get('query')
        sym = r.get('symbol')
        name = r.get('name')
        entrez = r.get('entrezgene')
        ensg = None
        ensembl_field = r.get('ensembl')
        if isinstance(ensembl_field, dict):
            ensg = ensembl_field.get('gene')
        elif isinstance(ensembl_field, list) and ensembl_field:
            # tomar el primero si hay varios
            first = ensembl_field[0]
            if isinstance(first, dict):
                ensg = first.get('gene')

        # UniProt puede ser str, dict (Swiss-Prot/TrEMBL) o lista
        uniprot_vals = []
        up = r.get('uniprot')
        if isinstance(up, str):
            uniprot_vals = [up]
        elif isinstance(up, dict):
            for k in ('Swiss-Prot', 'TrEMBL'):
                v = up.get(k)
                if v:
                    if isinstance(v, list):
                        uniprot_vals.extend([str(x) for x in v])
                    else:
                        uniprot_vals.append(str(v))
        elif isinstance(up, list):
            uniprot_vals.extend([str(x) for x in up])

        uniprot_unique = sorted(set([u for u in uniprot_vals if isinstance(u, str)]))
        rows.append({
            "input": q,
            "symbol_official": sym,
            "gene_name": name,
            "entrez": entrez,
            "ensembl_gene": ensg,
            "uniprot_ids": ";".join(uniprot_unique) if uniprot_unique else ""
        })

    df = pd.DataFrame(rows)
    logging.info("MyGene mapping completado.")
    return df


# ------------------------ Paso 2: Biopython/Entrez (resúmenes opcionales) ------------------------ #

def fetch_ncbi_summaries(symbols: List[str], email: str) -> str:
    """
    Recupera un resumen de NCBI Gene para cada símbolo humano (texto).
    Devuelve un bloque de texto para guardar.
    """
    if not email:
        return "No se proporcionó email para NCBI; se omiten resúmenes.\n"

    Entrez.email = email
    lines = ["# Resúmenes NCBI Gene (Homo sapiens)\n"]
    for sym in symbols:
        try:
            query = f"{sym}[Gene] AND Homo sapiens[Organism]"
            handle = Entrez.esearch(db="gene", term=query)
            rec = Entrez.read(handle)
            handle.close()
            ids = rec.get("IdList", [])
            if not ids:
                lines.append(f"## {sym}\nNo encontrado en NCBI Gene.\n")
                time.sleep(0.3)
                continue
            gid = ids[0]
            # El formato de 'gb' para gene es texto; también se podría usar docsum/summary.
            handle = Entrez.efetch(db="gene", id=gid, rettype="gene_table", retmode="text")
            txt = handle.read()
            handle.close()
            lines.append(f"## {sym}\n{txt}\n")
            time.sleep(0.3)  # cortesía NCBI
        except Exception as e:
            lines.append(f"## {sym}\nError: {e}\n")
    return "\n".join(lines)


# ------------------------ Paso 3: GOATOOLS (descargas, background y ORA) ------------------------ #

GO_OBO_URL = "http://current.geneontology.org/ontology/go-basic.obo"
GOA_GAF_GZ_URL = "http://current.geneontology.org/annotations/goa_human.gaf.gz"


def download_file(url: str, dest_path: str):
    if os.path.exists(dest_path):
        logging.info(f"Ya existe: {dest_path} (se omite descarga)")
        return
    logging.info(f"Descargando: {url}")
    r = requests.get(url, timeout=120)
    r.raise_for_status()
    with open(dest_path, "wb") as f:
        f.write(r.content)
    logging.info(f"Guardado: {dest_path}")


def ensure_go_resources(data_dir: str) -> Tuple[str, str]:
    """
    Descarga go-basic.obo y goa_human.gaf (descomprime .gz) en data_dir.
    Retorna rutas (obo_path, gaf_path)
    """
    ensure_dir(data_dir)
    obo_path = os.path.join(data_dir, "go-basic.obo")
    gaf_gz_path = os.path.join(data_dir, "goa_human.gaf.gz")
    gaf_path = os.path.join(data_dir, "goa_human.gaf")

    download_file(GO_OBO_URL, obo_path)

    if not os.path.exists(gaf_path):
        download_file(GOA_GAF_GZ_URL, gaf_gz_path)
        logging.info(f"Descomprimiendo {gaf_gz_path} -> {gaf_path}")
        with gzip.open(gaf_gz_path, "rb") as f_in, open(gaf_path, "wb") as f_out:
            shutil.copyfileobj(f_in, f_out)
    else:
        logging.info(f"Ya existe: {gaf_path}")

    return obo_path, gaf_path


def build_background_from_gaf(gaf_path: str) -> List[str]:
    """
    Construye el universo (background) como TODOS los identificadores de la columna 2 del GAF (UniProtKB).
    """
    bg = set()
    with open(gaf_path, "r", encoding="utf-8", errors="ignore") as f:
        for line in f:
            if line.startswith("!"):
                continue
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 3:
                continue
            uniprot_id = cols[1]  # columna 2: UniProtKB AC o ID
            if uniprot_id:
                bg.add(uniprot_id)
    bg_list = sorted(bg)
    logging.info(f"Universo (background) GAF: {len(bg_list)} identificadores UniProt.")
    return bg_list


def run_goatools_enrichment(study_uniprot: List[str], obo_path: str, gaf_path: str, out_path: str):
    """
    Ejecuta ORA con GOATOOLS (Fisher + BH FDR), propagate_counts=True, y guarda resultados.
    """
    if not study_uniprot:
        logging.warning("Lista de estudio vacía para GOATOOLS; se omite.")
        return

    # Ontología GO
    logging.info("Cargando ontología GO (go-basic.obo)...")
    godag = obo_parser.GODag(obo_file=obo_path)

    # Asociaciones UniProt->GO
    logging.info("Cargando asociaciones (GAF)...")
    gene2go = read_gaf(gaf_path)

    # Universo
    background = build_background_from_gaf(gaf_path)

    # Ejecutar estudio
    logging.info("Ejecutando GOEnrichmentStudy (propagate_counts=True, FDR BH)...")
    go_study = GOEnrichmentStudy(
        background,
        gene2go,
        godag,
        propagate_counts=True,
        alpha=0.05,
        methods=['fdr_bh']
    )

    results = go_study.run_study(study_uniprot)

    # Convertir a DataFrame
    rows = []
    for r in results:
        rows.append({
            "GO_ID": r.GO,
            "name": r.name,
            "namespace": r.NS,  # BP/MF/CC
            "enrichment": r.enrichment,  # e/p
            "p_uncorrected": r.p_uncorrected,
            "p_fdr_bh": getattr(r, "p_fdr_bh", None),
            "study_count": r.study_count,
            "study_n": r.ratio_in_study[1] if r.ratio_in_study else None,
            "pop_count": r.pop_count,
            "pop_n": r.ratio_in_pop[1] if r.ratio_in_pop else None
        })
    df = pd.DataFrame(rows)
    # Orden sugerido: por FDR y luego por p crudo
    if not df.empty:
        df.sort_values(by=["p_fdr_bh", "p_uncorrected"], inplace=True, na_position="last")
    save_df(df, out_path)


# ------------------------ Paso 4: Enrichr (API) ------------------------ #

ENRICHR_ADD_URL = "https://maayanlab.cloud/Enrichr/addList"
ENRICHR_ENRICH_URL = "https://maayanlab.cloud/Enrichr/enrich"

def enrichr_enrichment(symbols: List[str],
                       libraries: List[str],
                       outdir: str) -> None:
    """
    Sube la lista a Enrichr y ejecuta ORA en las bibliotecas indicadas.
    Guarda una tabla por biblioteca con p-value, Z-score, combined score y (si viene) FDR/overlap.
    """
    if not symbols:
        logging.warning("Lista vacía para Enrichr; se omite.")
        return

    gene_str = "\n".join([s.upper() for s in symbols])
    files = {
        'list': (None, gene_str),
        'description': (None, 'Functional analysis gene list')
    }
    logging.info("Subiendo lista a Enrichr...")
    r = requests.post(ENRICHR_ADD_URL, files=files, timeout=60)
    if r.status_code != 200:
        logging.error(f"Error al subir lista a Enrichr: {r.status_code} {r.text}")
        return
    user_list_id = r.json().get("userListId")
    logging.info(f"Enrichr userListId={user_list_id}")

    for lib in libraries:
        params = {"userListId": user_list_id, "backgroundType": lib}
        logging.info(f"Consultando Enrichr: {lib}")
        r = requests.get(ENRICHR_ENRICH_URL, params=params, timeout=60)
        if r.status_code != 200:
            logging.error(f"Error en Enrichr ({lib}): {r.status_code} {r.text}")
            continue
        data = r.json()
        # data[lib] es una lista de filas (lista posicional)
        # Índices más habituales (según docs de Enrichr):
        # 0: rank, 1: term_name, 2: pvalue, 3: zscore, 4: combined_score,
        # 5: overlapping_genes, 6: adjusted_pvalue (FDR), etc. (puede variar)
        rows = []
        for rec in data.get(lib, []):
            term_name = rec[1] if len(rec) > 1 else ""
            pvalue = rec[2] if len(rec) > 2 else None
            zscore = rec[3] if len(rec) > 3 else None
            combined = rec[4] if len(rec) > 4 else None
            overlap = rec[5] if len(rec) > 5 else ""
            adjp = rec[6] if len(rec) > 6 else None
            rows.append({
                "library": lib,
                "term": term_name,
                "pvalue": pvalue,
                "p_adj_FDR": adjp,
                "z_score": zscore,
                "combined_score": combined,
                "overlap_genes": overlap
            })
        df = pd.DataFrame(rows)
        out_path = os.path.join(outdir, f"enrichr_{lib}.tsv")
        if not df.empty:
            df.sort_values(by=["p_adj_FDR", "pvalue", "combined_score"], inplace=True, na_position="last")
        save_df(df, out_path)


# ------------------------ Paso 5: STRING (API) ------------------------ #

STRING_API_BASE = "https://version-11-5.string-db.org/api"
def string_enrichment(symbols: List[str], species: int, caller_identity: str, out_path: str):
    """
    Ejecuta el endpoint 'enrichment' de STRING v11.5 con una lista de símbolos humanos.
    Guarda todas las categorías devueltas, incluyendo FDR y lista de genes por término.
    """
    if not symbols:
        logging.warning("Lista vacía para STRING; se omite.")
        return

    method = "enrichment"
    url = "/".join([STRING_API_BASE, "json", method])
    params = {
        "identifiers": "%0d".join(symbols),
        "species": species,
        "caller_identity": caller_identity or "functional_analysis_script"
    }
    logging.info("Consultando STRING (enrichment)...")
    r = requests.post(url, data=params, timeout=120)
    if r.status_code != 200:
        logging.error(f"Error STRING enrichment: {r.status_code} {r.text}")
        return
    try:
        data = r.json()
    except Exception:
        data = json.loads(r.text)

    rows = []
    for row in data:
        rows.append({
            "term": row.get("term"),
            "category": row.get("category"),
            "description": row.get("description"),
            "fdr": row.get("fdr"),
            "p_value": row.get("p_value"),
            "preferredNames": ",".join(row.get("preferredNames", [])),
            "inputGenes": ",".join(row.get("inputGenes", [])),
        })
    df = pd.DataFrame(rows)
    if not df.empty:
        df.sort_values(by=["fdr", "p_value"], inplace=True, na_position="last")
    save_df(df, out_path)


# ------------------------ Empaquetado CLI ------------------------ #

def parse_args():
    p = argparse.ArgumentParser(
        description="Tarea 1: Análisis funcional (COX4I2, ND1, ATP6) con MyGene, GOATOOLS, Enrichr y STRING."
    )
    p.add_argument("--input", required=True, help="Ruta del archivo de genes de entrada (uno por línea).")
    p.add_argument("--outdir", required=True, help="Directorio de salida para resultados.")
    p.add_argument("--data_dir", default="data", help="Directorio para recursos GO (obo/gaf).")
    p.add_argument("--species", default="human", help="Especie para MyGene (por defecto: human).")
    p.add_argument("--taxid", type=int, default=9606, help="NCBI taxid para STRING (por defecto 9606).")
    p.add_argument("--email_topt", default="", help="Email para Biopython/Entrez (opcional).")
    p.add_argument("--string_id", default="", help="Caller identity para STRING (email/ID).")
    p.add_argument("--enrichr_libs", default="GO_Biological_Process_2023,Reactome_2022,KEGG_2019_Human",
                   help="Lista de bibliotecas Enrichr separadas por coma.")
    return p.parse_args()


def main():
    setup_logging()
    args = parse_args()

    ensure_dir(args.outdir)

    # 1) Leer genes de entrada
    genes_input = read_gene_list(args.input)
    logging.info(f"Genes de entrada: {genes_input}")

    # 2) MyGene: mapping
    df_map = mygene_mapping(genes_input, species=args.species)
    save_df(df_map, os.path.join(args.outdir, "id_mapping.tsv"))

    # Obtener símbolos oficiales (para Enrichr/STRING)
    symbols_official = [s for s in df_map["symbol_official"].dropna().astype(str).tolist() if s]
    # Para GOATOOLS: UniProt IDs
    uniprot_pool = []
    for ids in df_map["uniprot_ids"].fillna("").tolist():
        if ids:
            uniprot_pool.extend([x.strip() for x in ids.split(";") if x.strip()])
    study_uniprot = sorted(set(uniprot_pool))
    logging.info(f"Símbolos oficiales (para Enrichr/STRING): {symbols_official}")
    logging.info(f"UniProt para GOATOOLS: {study_uniprot}")

    # 3) (Opcional) Resúmenes NCBI Gene
    if args.email_topt:
        summaries = fetch_ncbi_summaries(symbols_official, args.email_topt)
        summaries_path = os.path.join(args.outdir, "ncbi_gene_summaries.txt")
        with open(summaries_path, "w", encoding="utf-8") as f:
            f.write(summaries)
        logging.info(f"Guardado: {summaries_path}")
    else:
        logging.info("Sin email para Entrez; se omiten resúmenes NCBI.")

    # 4) GOATOOLS: descargar recursos y ejecutar ORA
    obo_path, gaf_path = ensure_go_resources(args.data_dir)
    goatools_out = os.path.join(args.outdir, "goatools_enrichment.tsv")
    run_goatools_enrichment(study_uniprot, obo_path, gaf_path, goatools_out)

    # 5) Enrichr (ORA)
    enrichr_libs = [x.strip() for x in args.enrichr_libs.split(",") if x.strip()]
    enrichr_enrichment(symbols_official, enrichr_libs, args.outdir)

    # 6) STRING (ORA)
    string_out = os.path.join(args.outdir, "string_enrichment.tsv")
    string_enrichment(symbols_official, args.taxid, args.string_id, string_out)

    logging.info("Análisis funcional completado ✅")


if __name__ == "__main__":
    try:
        main()
    except requests.exceptions.RequestException as e:
            logging.error(f"Error de red: {e}")
            sys.exit(1)
    except Exception as e:
        logging.error(f"Error inesperado: {e}")
        sys.exit(1)

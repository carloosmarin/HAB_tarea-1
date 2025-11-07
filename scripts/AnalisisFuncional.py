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
- STRINGdb (API) para ORA con respaldo de red de PPIs de STRING.

Flujo:
1) Leer genes de entrada (símbolos o IDs).
2) Normalizar/validar y mapear a UniProt/Ensembl/Entrez con mygene.
3) (Opcional) Descargar resúmenes de NCBI Gene con Biopython/Entrez.
4) Descargar go-basic.obo y goa_human.gaf; construir universo (background) y asociaciones (UniProt→GO).
5) ORA de GO con GOATOOLS (propagate_counts=True; FDR BH).
6) ORA con Enrichr (API) en bibliotecas seleccionadas (por defecto: GO BP, Reactome, KEGG).
7) ORA con STRING (API) y guardar resultados.
8) Guardar outputs en 'results/' y limpieza de archivos pesados.

CLI (ejemplo PowerShell):
    python scripts/AnalisisFuncional.py `
      --input data/genes_input.txt `
      --outdir results `
      --email_topt TU_EMAIL@DOMINIO `
      --string_id TU_EMAIL@DOMINIO `
      --input_id_type auto `
      --go_bg_size 2000

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
import re
import random
from time import sleep
from typing import List, Tuple

import requests
import pandas as pd
import mygene
from Bio import Entrez  # Biopython (opcional)

# GOATOOLS
from goatools import obo_parser
from goatools.associations import read_gaf
from goatools.go_enrichment import GOEnrichmentStudy


# ------------------------ Helpers de logging / sistema ------------------------ #

def setup_logging():
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s [%(levelname)s] %(message)s",
        datefmt="%H:%M:%S"
    )


def ensure_dir(path: str):
    os.makedirs(path, exist_ok=True)


# ------------------------ Detección / normalización de IDs ------------------------ #

ENSEMBL_RE = re.compile(r"^ENSG\d{11}$", re.IGNORECASE)
UNIPROT_RE = re.compile(r"^[A-NR-Z0-9]{6,10}$")  # heurística simple
ENTREZ_RE = re.compile(r"^\d+$")


def detect_id_type(s: str) -> str:
    s = s.strip()
    if ENSEMBL_RE.match(s):
        return "ensembl"
    if ENTREZ_RE.match(s):
        return "entrez"
    # Si parece uniprot y está en mayúsculas
    if UNIPROT_RE.match(s) and s.upper() == s:
        return "uniprot"
    return "symbol"


def fix_human_mito_aliases(symbols: List[str]) -> List[str]:
    """
    Corrige alias frecuentes de genes mitocondriales humanos a los símbolos HGNC oficiales.
    Ej.: ND1 -> MT-ND1, ATP6 -> MT-ATP6, etc.
    """
    mito_map = {
        "ND1": "MT-ND1",
        "ND2": "MT-ND2",
        "ND3": "MT-ND3",
        "ND4": "MT-ND4",
        "ND4L": "MT-ND4L",
        "ND5": "MT-ND5",
        "ND6": "MT-ND6",
        "ATP6": "MT-ATP6",
        "ATP8": "MT-ATP8",
        "COX1": "MT-CO1",
        "COX2": "MT-CO2",
        "COX3": "MT-CO3",
        "CYTB": "MT-CYB",
    }
    fixed = []
    for s in symbols:
        s_clean = s.strip().upper().replace(",", "")
        fixed.append(mito_map.get(s_clean, s_clean))
    return fixed


def read_gene_list(input_path: str) -> List[str]:
    """
    Lee el archivo de genes de entrada. Admite formato con comas (COX4I2, ND1, ATP6)
    o formato por líneas (uno por línea).
    """
    genes = []
    with open(input_path, "r", encoding="utf-8") as f:
        text = f.read().strip()
        text = text.replace(",", "\n")
        for line in text.splitlines():
            g = line.strip()
            if g and not g.startswith("#"):
                genes.append(g)
    if not genes:
        raise ValueError(f"El archivo de entrada '{input_path}' no contiene genes válidos.")
    return genes


def save_df(df: pd.DataFrame, path: str):
    if df is None or df.empty:
        logging.warning(f"No hay datos para guardar en {path}")
        return
    df.to_csv(path, sep="\t", index=False)
    logging.info(f"Guardado: {path}")


# ------------------------ Red robusta (reintentos/backoff) ------------------------ #

def robust_request(method, url, *, max_retries=3, backoff=2.0, timeout=60, **kwargs):
    """
    Envuelve requests.get/post con reintentos exponenciales y logging claro.
    method: 'get' o 'post'
    """
    func = requests.get if method.lower() == "get" else requests.post
    last_err = None
    for attempt in range(1, max_retries + 1):
        try:
            resp = func(url, timeout=timeout, **kwargs)
            resp.raise_for_status()
            return resp
        except requests.exceptions.RequestException as e:
            last_err = e
            logging.warning(f"[NET] {method.upper()} {url} falló (intento {attempt}/{max_retries}): {e}")
            if attempt < max_retries:
                sleep(backoff ** attempt)
    raise last_err


# ------------------------ Paso 1: MyGene (normalización / mapping) ------------------------ #

def mygene_mapping(genes: List[str],
                   species: str = "human",
                   input_id_type: str = "auto",
                   output_cols: List[str] | None = None) -> Tuple[pd.DataFrame, List[str]]:
    """
    Normaliza y convierte IDs con MyGene:
    - Admite entrada: symbol|entrez|ensembl|uniprot|auto
    - Devuelve DataFrame con columnas de salida y lista de no mapeados.
    """
    if output_cols is None:
        output_cols = ["symbol_official", "gene_name", "entrez", "ensembl_gene", "uniprot_ids"]

    scopes_map = {
        "symbol": "symbol,alias",
        "entrez": "entrezgene",
        "ensembl": "ensembl.gene",
        "uniprot": "uniprot",
    }

    if input_id_type == "auto":
        detected = {detect_id_type(g) for g in genes}
        if len(detected) == 1:
            scopes = scopes_map[list(detected)[0]]
        else:
            scopes = "symbol,alias,entrezgene,ensembl.gene,uniprot"
    else:
        scopes = scopes_map.get(input_id_type, "symbol,alias")

    mg = mygene.MyGeneInfo()
    logging.info(f"Consultando MyGene para {len(genes)} genes... (scopes='{scopes}')")
    try:
        res = mg.querymany(
            genes,
            scopes=scopes,
            fields="symbol,name,entrezgene,ensembl.gene,uniprot",
            species=species,
            as_dataframe=False
        )
    except Exception as e:
        raise RuntimeError(f"Fallo consultando MyGene: {e}")

    rows, missing = [], []
    for r in res:
        if r.get('notfound'):
            missing.append(r.get('query'))
            continue

        q = r.get('query')
        sym = r.get('symbol')
        name = r.get('name')
        entrez = r.get('entrezgene')

        ensg = None
        ensembl_field = r.get('ensembl')
        if isinstance(ensembl_field, dict):
            ensg = ensembl_field.get('gene')
        elif isinstance(ensembl_field, list) and ensembl_field:
            first = ensembl_field[0]
            if isinstance(first, dict):
                ensg = first.get('gene')

        # UniProt puede ser str, dict o lista
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
    # Asegura columnas solicitadas
    for c in ["input"] + output_cols:
        if c not in df.columns:
            df[c] = ""

    if missing:
        logging.warning(f"MyGene: {len(missing)} no mapeados: {missing}")
    else:
        logging.info("MyGene mapping completado (sin faltantes).")

    return df, missing


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
    r = robust_request("get", url, timeout=120)
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


def run_goatools_enrichment(study_uniprot: List[str], obo_path: str, gaf_path: str, out_path: str, bg_size: int):
    """
    Ejecuta ORA con GOATOOLS (Fisher + BH FDR), propagate_counts=True, y guarda resultados significativos.
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

    # Universo (background) reducido
    background_full = build_background_from_gaf(gaf_path)

    random.seed(42)
    background = random.sample(background_full, min(bg_size, len(background_full)))

    # Aseguramos que los genes de estudio estén incluidos
    for gene in study_uniprot:
        if gene not in background:
            background.append(gene)

    logging.info(f"Usando background reducido de {len(background)} genes (de {len(background_full)} totales).")

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

    # Filtrar resultados significativos (FDR < 0.05)
    significant = [r for r in results if getattr(r, "p_fdr_bh", 1.0) < 0.05]

    if not significant:
        logging.info("No se encontraron términos GO significativos.")
        # No guardamos archivo vacío; el usuario verá el warning de save_df si se quisiera escribir.
        return

    logging.info(f"Se encontraron {len(significant)} términos GO significativos (FDR < 0.05).")

    # Convertir los resultados significativos a DataFrame compacto
    rows = []
    for r in significant:
        rows.append({
            "GO_ID": r.GO,
            "name": r.name,
            "namespace": r.NS,
            "enrichment": r.enrichment,
            "p_uncorrected": r.p_uncorrected,
            "p_fdr_bh": getattr(r, "p_fdr_bh", None),
            "study_count": r.study_count,
            "study_n": r.ratio_in_study[1] if r.ratio_in_study else None,
            "pop_count": r.pop_count,
            "pop_n": r.ratio_in_pop[1] if r.ratio_in_pop else None
        })
    df = pd.DataFrame(rows)
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
    r = robust_request("post", ENRICHR_ADD_URL, files=files, timeout=60)
    user_list_id = r.json().get("userListId")
    logging.info(f"Enrichr userListId={user_list_id}")

    for lib in libraries:
        params = {"userListId": user_list_id, "backgroundType": lib}
        logging.info(f"Consultando Enrichr: {lib}")
        r = robust_request("get", ENRICHR_ENRICH_URL, params=params, timeout=60)
        data = r.json()
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
    r = robust_request("post", url, data=params, timeout=120)

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
    p.add_argument("--input", required=True, help="Ruta del archivo de genes de entrada (con comas o por líneas).")
    p.add_argument("--outdir", required=True, help="Directorio de salida para resultados.")
    p.add_argument("--data_dir", default="data", help="Directorio para recursos GO (obo/gaf).")
    p.add_argument("--species", default="human", help="Especie para MyGene (por defecto: human).")
    p.add_argument("--taxid", type=int, default=9606, help="NCBI taxid para STRING (por defecto 9606).")
    p.add_argument("--email_topt", default="", help="Email para Biopython/Entrez (opcional).")
    p.add_argument("--string_id", default="", help="Caller identity para STRING (email/ID).")
    p.add_argument("--enrichr_libs", default="GO_Biological_Process_2023,Reactome_2022,KEGG_2019_Human",
                   help="Lista de bibliotecas Enrichr separadas por coma.")
    p.add_argument("--input_id_type", choices=["auto", "symbol", "entrez", "ensembl", "uniprot"],
                   default="auto",
                   help="Tipo de ID de entrada (auto|symbol|entrez|ensembl|uniprot).")
    p.add_argument("--output_id_cols", default="symbol_official,entrez,ensembl_gene,uniprot_ids",
                   help="Columnas de ID a asegurar en el mapeo (coma).")
    p.add_argument("--go_bg_size", type=int, default=2000,
                   help="Tamaño del background aleatorio para GOATOOLS (por defecto 2000).")
    return p.parse_args()


def main():
    setup_logging()
    args = parse_args()

    ensure_dir(args.outdir)

    # 1) Leer genes de entrada
    genes_input = read_gene_list(args.input)
    logging.info(f"Genes de entrada: {genes_input}")
    # Normalizar alias mitocondriales (humano) si el input es simbólico/auto
    if args.input_id_type in ("auto", "symbol"):
        genes_input = fix_human_mito_aliases(genes_input)
    logging.info(f"Genes normalizados: {genes_input}")

    # 2) MyGene: mapping con detección/selección de ID de entrada
    output_cols = [c.strip() for c in args.output_id_cols.split(",") if c.strip()]
    df_map, missing = mygene_mapping(
        genes_input,
        species=args.species,
        input_id_type=args.input_id_type,
        output_cols=output_cols
    )
    save_df(df_map, os.path.join(args.outdir, "id_mapping.tsv"))

    if missing:
        missing_path = os.path.join(args.outdir, "id_unmapped.txt")
        with open(missing_path, "w", encoding="utf-8") as fh:
            fh.write("\n".join(missing) + "\n")
        logging.warning(f"Guardado reporte de no mapeados: {missing_path}")

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

    if not symbols_official:
        logging.error("No se obtuvieron símbolos oficiales tras el mapeo. Revisa el tipo de ID de entrada o los nombres.")
        sys.exit(2)

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
    run_goatools_enrichment(study_uniprot, obo_path, gaf_path, goatools_out, bg_size=args.go_bg_size)

    # 5) Enrichr (ORA)
    enrichr_libs = [x.strip() for x in args.enrichr_libs.split(",") if x.strip()]
    enrichr_enrichment(symbols_official, enrichr_libs, args.outdir)

    # 6) STRING (ORA)
    string_out = os.path.join(args.outdir, "string_enrichment.tsv")
    string_enrichment(symbols_official, args.taxid, args.string_id, string_out)

    # 7) Limpieza: eliminar archivos pesados intermedios
    try:
        for file_to_remove in [
            os.path.join(args.data_dir, "goa_human.gaf"),
            os.path.join(args.data_dir, "goa_human.gaf.gz")
        ]:
            if os.path.exists(file_to_remove):
                os.remove(file_to_remove)
                logging.info(f"Archivo eliminado para optimizar espacio: {file_to_remove}")
    except Exception as e:
        logging.warning(f"No se pudo eliminar archivo intermedio: {e}")

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

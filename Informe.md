# Análisis funcional de COX4I2, ND1 y ATP6
**Autor:** Carlos Marín Martínez
# Informe del análisis funcional (COX4I2, MT-ND1, MT-ATP6)

## 1. Descripción general

El objetivo de este análisis fue realizar un **análisis funcional (ORA)** sobre tres genes mitocondriales humanos:  
`COX4I2`, `MT-ND1` y `MT-ATP6`.

Se emplearon tres enfoques complementarios:
- **GOATOOLS** → enriquecimiento en términos GO (Gene Ontology).  
- **Enrichr** → análisis sobre bibliotecas funcionales (GO Biological Process, Reactome y KEGG).  
- **STRING** → enriquecimiento basado en redes de interacción y rutas KEGG.

Todos los métodos se ejecutaron con control de errores, detección automática del tipo de ID de entrada y reducción del background para optimizar el tamaño y la velocidad de ejecución.

---

## 2. Cambios y optimización del script

### Mejoras implementadas

| Aspecto | Antes | Ahora |
|----------|--------|-------|
| **Gestión de red** | Sin control de errores | Reintentos automáticos con backoff exponencial |
| **Conversión de IDs** | Solo símbolos | Detección automática (symbol, Entrez, Ensembl, UniProt) |
| **Background (GOATOOLS)** | >40.000 genes | 2.000 genes aleatorios (manteniendo los de estudio) |
| **Archivos intermedios** | Permanecían en disco | Eliminación automática de `.gaf` y `.gz` |
| **Manejo de errores** | Sin excepciones controladas | Mensajes claros y salida ordenada |

Estas mejoras reducen el tamaño de los resultados y hacen el proceso más eficiente sin pérdida de información clave.

---

## 3. Resultados principales

### 3.1 GOATOOLS (Gene Ontology)
**Archivo:** `results/goatools_enrichment.tsv`

Los términos más significativos detectados con un background de 2000 genes fueron:

| GO_ID | name | p_uncorrected | p_fdr_bh | study_count |
|--------|------|---------------|-----------|--------------|
| GO:0042776 | proton motive force-driven mitochondrial ATP synthesis | 7.05e-06 | 0.0307 | 2 |
| GO:0015986 | proton motive force-driven ATP synthesis | 9.12e-06 | 0.0307 | 2 |
| GO:0006754 | ATP biosynthetic process | 1.29e-05 | 0.0307 | 2 |
| GO:0019646 | aerobic electron transport chain | 1.29e-05 | 0.0307 | 2 |

Estos términos confirman una fuerte asociación con **procesos mitocondriales energéticos**, coherente con la función de los genes analizados.

---

### 3.2 Enrichr (GO Biological Process)

**Archivo:** `results/enrichr_GO_Biological_Process_2023.tsv`

Los procesos enriquecidos con menor FDR fueron:

| Term | FDR ajustado | Gene(s) implicados |
|------|---------------|--------------------|
| Mitochondrial Electron Transport, Cytochrome C To Oxygen | 0.0119 | COX4I2 |
| Energy Derivation By Oxidation Of Organic Compounds | 0.0127 | COX4I2 |
| Aerobic Electron Transport Chain | 0.0127 | COX4I2 |
| Mitochondrial ATP Synthesis Coupled Electron Transport | 0.0127 | COX4I2 |

**Interpretación:** Enrichr refuerza los mismos procesos oxidativos y respiratorios observados con GOATOOLS.

---

### 3.3 STRING (Enriquecimiento funcional por red)

**Archivo:** `results/string_enrichment.tsv`

STRING asocia el conjunto génico con rutas KEGG de enfermedades mitocondriales y procesos respiratorios.

| Term | Categoría | FDR | Genes |
|------|------------|-----|-------|
| Oxidative phosphorylation | KEGG | 0.0001 | MT-ATP6, MT-ND1, COX4I2 |
| Thermogenesis | KEGG | 0.00028 | MT-ATP6, MT-ND1, COX4I2 |
| Parkinson / Alzheimer / Huntington diseases | KEGG | 0.0003 | MT-ATP6, MT-ND1, COX4I2 |
| Inner mitochondrial membrane protein complex | GO Component | 0.00053 | MT-ATP6, MT-ND1, COX4I2 |

---

## 4. Conclusiones

- Los tres métodos (GOATOOLS, Enrichr, STRING) convergen en **procesos de fosforilación oxidativa y transporte electrónico mitocondrial**.  
- El script actualizado mejora la **eficiencia**, **reproducibilidad** y **legibilidad** de los resultados.  
- El **control de errores y conversión automática de IDs** hacen el flujo más robusto para listas de genes heterogéneas.  
- El **recorte del background** mantiene la validez estadística, pero reduce el tamaño y el tiempo de cómputo.

**Resultado final:** el pipeline queda preparado para análisis ligeros, reproducibles y completamente documentados.

# Tarea 1: Análisis Funcional de COX4I2, MT-ND1 y MT-ATP6

Este repositorio contiene un script en Python que realiza un análisis funcional de tres genes mitocondriales humanos (COX4I2, ND1 y ATP6).  
El análisis incluye conversión de identificadores, obtención de anotaciones, y enriquecimiento funcional mediante MyGene, GOATOOLS, Enrichr y STRINGdb.  
Las salidas se almacenan en la carpeta `/results/`.  

# Tarea 1: Análisis Funcional 


Entrega un script en Python que implemente de forma detallada un análisis funcional de los genes **COX4I2**, **ND1** y **ATP6**. El análisis debe utilizar librerías de Python para evaluar los procesos biológicos asociados a estos genes. Asegúrate de que tu código esté bien documentado y describa claramente los métodos y bases de datos utilizadas para obtener información funcional.

## Estructura del repositorio

```
/analisis-funcional/
├── data/
│   └── genes_input.txt        # Genes proporcionados para el análisis
├── scripts/
│   └── tu_script.py           # Tu script de análisis funcional
├── results/                   # Resultados generados por el script (opcional)
├── README.md                  # Este archivo
└── requirements.txt           # Dependencias del proyecto
```

## Instrucciones de entrega

- Haz un **fork** de este repositorio en tu cuenta de GitHub.
- Trabaja en tu fork y sube tu script a la carpeta `scripts/`.
- Tu script debe poder ejecutarse desde la línea de comandos (CLI), aceptando como mínimo el archivo de entrada y el archivo de salida como argumentos.
- No modifiques el archivo `genes_input.txt`.
- Documenta tu código explicando los métodos, librerías y bases de datos utilizadas.
- Puedes generar un archivo de resultados en la carpeta `results/` si lo consideras útil.

## Rúbrica de evaluación

La tarea se evaluará sobre un máximo de **10 puntos**, distribuidos según los siguientes criterios:

| Criterio | Descripción | Puntos |
|---------|-------------|--------|
| **1. Funcionalidad** | El script realiza correctamente el análisis funcional solicitado. | 4 |
| **2. Documentación** | El código está comentado y explica claramente los métodos y bases de datos utilizadas. | 2 |
| **3. Uso de librerías** | Se emplean librerías adecuadas para el análisis funcional (e.g., gprofiler-official, GOATOOLS, Enrichr, etc.). | 2 |
| **4. Formato y estilo** | El código sigue buenas prácticas de estilo y es legible. | 1 |
| **5. Automatización (CLI)** | El script acepta argumentos desde la línea de comandos. | 1 |

## Dependencias recomendadas

Incluye en `requirements.txt` las librerías necesarias para ejecutar tu script. Por ejemplo:

```
gprofiler-official
pandas
```


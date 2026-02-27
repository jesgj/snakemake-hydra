#!/bin/bash
#SBATCH --job-name=skm # Nombre del trabajo
#SBATCH --output=skm-%j.out # Archivo de salida
#SBATCH --error=skm-%j.err # Archivo de errores
#SBATCH --partition=standard # Partición donde ejecutar el trabajo
#SBATCH --nodes=1 # Número de nodos
#SBATCH --ntasks-per-node=1 # Número de procesos/tareas por nodo
#SBATCH --cpus-per-task=22 # Número de núcleos por proceso/tarea
#SBATCH --mem=26G
#SBATCH --time=06:00:00 # Tiempo máximo de ejecución (horas:minutos:segundos)


wdir=/home/jlgarcia/CTA_TFM_UOC
cd $wdir
pwd

# run snakemake
pixi run snakemake --cores 22 -s Snakefile --until chip_cr_fastqc_raw_pe --latency-wait 60 --printshellcmds

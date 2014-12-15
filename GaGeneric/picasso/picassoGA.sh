#!/bin/bash
# MASTER
# numero de cpus que empleara el calculo:
#PBS -l ncpus=1
# memoria que empleara el calculo:
#PBS -l mem=400mb
# como mucho tardara X horas
#PBS -l walltime=10:00:00
# si se quiere usar un array job, ahora mismo esta comentado:
# nada PBS -J 1-1
# para que lo envie al superdome:
#PBS -q routex86

# para que vaya al directorio actual:
cd "$PBS_O_WORKDIR"

# programa a ejecutar, con sus argumentos:
matlab -nodisplay -nojvm -r picasso1 > picasso1Output.txt


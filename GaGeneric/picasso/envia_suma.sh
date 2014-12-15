#!/bin/bash

# numero de cpus que empleara el calculo:
#PBS -l ncpus=1
# memoria que empleara el calculo:
#PBS -l mem=2gb
# como mucho tardara 5 horas
#PBS -l walltime=5:00:00
# si se quiere usar un array job, ahora mismo esta comentado:
#PBS -J 1-4
# para que lo envie al superdome:
#PBS -q routex86

# para que vaya al directorio actual:
cd "$PBS_O_WORKDIR"

# programa a ejecutar, con sus argumentos:


echo 2+$PBS_ARRAY_INDEX > suma_${PBS_ARRAY_INDEX}.m

# programa a ejecutar, con sus argumentos:

time matlab < suma_${PBS_ARRAY_INDEX}.m


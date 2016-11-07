#!/bin/bash
#
#$ -cwd
#$ -j y
#$ -S /bin/bash
#$ -pe openmpi 5
#$ -V
#

     mpirun -n 5 ./nhwave

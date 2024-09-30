#!  /usr/bin/bash

. /local/env/envnextflow-23.10.0.sh

nextflow run Contams.nf -profile slurm,conda
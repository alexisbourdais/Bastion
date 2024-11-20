#!  /usr/bin/bash

### Genouest

nextflow run Bastion.nf \
-profile standard,singularity \
--workflow analysis \
--singularity "-B /home:/home"

#--workflow analysis \
#--workflow setup \

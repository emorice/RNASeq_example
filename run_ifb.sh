#!/bin/bash
DIR='data/mydatalocal'
rsync -tu snakefile ensure_deps.sh ubuntu@134.158.247.94:$DIR
ssh ubuntu@134.158.247.94 bash -lc \'cd $DIR \; conda activate \; ./ensure_deps.sh \; snakemake --rerun-incomplete -j all\'

#!/bin/bash
DIR='data/mydatalocal/EMT'
rsync -tu snakefile ensure_deps.sh ubuntu@134.158.247.94:$DIR
if [ "$1" = "unlock" ]; then
	ssh ubuntu@134.158.247.94 bash -lc \'cd $DIR \|\| \{ mkdir $DIR \; cd $DIR \;\}\; conda activate \; ./ensure_deps.sh \; snakemake --unlock\'
fi
ssh ubuntu@134.158.247.94 bash -lc \'cd $DIR \|\| \{ mkdir $DIR \; cd $DIR \;\}\; conda activate \; ./ensure_deps.sh \; snakemake --rerun-incomplete -j all\'

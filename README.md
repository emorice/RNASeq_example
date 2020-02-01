# NGS example pipelines

Each directory contains a pipeline divided into three files:
 * a `snakefile` containing the pipeline, that can be run with `snakemake` or `snakemake -j <n>` for `n` cores
 * a `ensure_deps.sh` script that checks and if neccessary installs the required tools from conda, including snakemake. It is designed to be reentrant, if one of the tools is already available, no install will be attempted.
 * a `run_ifb.sh` script demonstrating how to run the pipeline from scratch on a distant machine. It takes cares of copying the pipeline script, installing dependencies, creating working directory and running the pipeline. The provided host will of course not work for you and is kept as an example.
 
 ## Input files
 Each pipeline requires input files :
  * For EMT, a `source` subdirectory containing the `.fastq` files, and a set of `chr18.fa`, `chr18.gtf`, `chr18.idx/`annotated and STAR-indexed genome files.
  * For Variants, an `exome` directory containing the `.fastq.gz` file and a set of `chr16.fa.gz`, `chr16.fa.gz.*`, `gencode.v24lift37.basic.annotation.gtf.gz` annotated and BWA-indexed genome files. Warning: `chr16.fa.gz` must be recompressed from gzip to BGZF if necessary.
  
 For convenience, the pipeline are written to fetch, extract and index these files autonomously, but due to the incremental and reentrant nature of snakemake pipelines, these steps will **not** be run if the files are already present, meaning that you can prepopulate the working directory with your input files, following the naming conventions, and the input steps will be automatically skipped.
 
 ## Output
 The EMT pipeline produces a `counts.txt` file with the read count per gene.
 The Variants pipeline produces a `variants.features` and `variants.*.vcf` files with the called and bed-annotated variants.
 Fastqc steps are also included, they are run automatically for the Variant pipeline and results stored in `quality/`.
 
 ## Snakemake assets
 Compared to a bash script, writting these pipelines with snakemake:
  * Eases development and extension: one can simply add new rules (steps) in the script, change the objective rule, and re-run: the unchanged steps will be skipped to go straight to the new steps, while still producing a pipeline able to run form start to finish on another machine.
  * Improves modularity: each step has explicit inputs and outputs, making it easier to reuse and adapt them to another pipeline.
  * Improves parallelism: independant steps will be run in parallel. For instance, in the Variant pipeline, both normal and tumor pileups will run in parallel, using two cores, in the EMT pipeline, a file can be run both trimmed and run through fastqc at the same time.
  * Improves scaling: snakemake natively supports remote execution on a cluster, so one can for example run `snakemake --cluster qsub -j 32` to run the dteps in a distributed fashion on a SGE platform, for instance, without changing the pipeline script.
  

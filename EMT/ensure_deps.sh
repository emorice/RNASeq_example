#!/bin/bash
{ fastqc -v &> /dev/null ;} || conda install -y fastqc
{ trimmomatic -version &> /dev/null ;} || conda install -y trimmomatic
{ STAR --version &> /dev/null ;} || conda install -y STAR
{ snakemake -v &> /dev/null ;} || conda install -y snakemake
{ samtools help &> /dev/null ;} || conda install -y samtools
{ featureCounts -v &> /dev/null ;} || conda install -y subread

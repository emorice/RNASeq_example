#!/bin/bash
{ fastqc -v &> /dev/null ;} || conda install -y fastqc
{ trimmomatic -version &> /dev/null ;} || conda install -y trimmomatic
{ bwa |& grep 'Program: bwa' &> /dev/null ;} || conda install -y bwa
{ snakemake -v &> /dev/null ;} || conda install -y snakemake
{ samtools help &> /dev/null ;} || conda install -y samtools
{ varscan -v &> /dev/null ;} || conda install -y varscan
{ bedtools --version &> /dev/null ;} || conda install -y bedtools


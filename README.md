# nf_ccle

Sam Lee
email: lee.sa@wehi.edu.au


## about

some info about the data: https://www.ebi.ac.uk/biostudies/files/E-MTAB-2770/E-MTAB-2770.idf.txt




## To run

```bash
# start a tmux / screen session

module load nextflow

nextflow run main.nf -c ccle_params.config -resume -with-dag flowchart.mmd
```


## Analysis steps

Using (https://github.com/samleenz/nf_modules) for analysis code

- Quality control
  - fastqc
  - multiQC
- Read Processing
  - Trimming with cutadapt
- Alignment
  - STAR two-pass alignment
- Quantification  
  - featureCounts 

## Versions

- Ensembl release 108 for txome
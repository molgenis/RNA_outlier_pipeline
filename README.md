# RNA outlier pipeline
Used in manuscript "Accelerating rare disease diagnostics by explainable linking of DNA and RNA using a simplified and interactive RNA-guided workflow"

This repository adheres to the following structure

```
├───Nextflow
│       ├───MAE
│       ├───fraser
│       ├───html_report
│       ├───outrider
│       ├───main.nf
│       ├───nextflow.config
├───LICENSE
└───README.md
```

## Nextflow pipeline

In the Nextflow directory, you'll find the code used to create a RNA outlier pipeline (nextflow version 24.04.2) to generate a .tsv files containing aberrantly expressed genes and a .tsv containing aberrantly spliced genes. This pipeline is based on [OUTRIDER](version 1.20.1) and [FRASER](version 1.99.4) (see "References")

### Requirements

- GNU-based Linux (Rocky Linux 9.3 used in manuscript)
- Bash ≥ 3.2
- 8 CPUs
- 64GB RAM
- [Nextflow(version 24.04.4)](https://www.nextflow.io/)
- [Apptainer(version 1.3.0-1.el9)](https://apptainer.org/)
- Python ≥ 3.10.1
  - Packages:
    - json
    - argparse
    - pandas
    - plotly.express
    - numpy
    - jinja2
- Singularity container with OUTRIDER and FRASER
  - Add the Sylabscloud to the remotes on Apptainer if not already configured
  <code>
  apptainer remote add --no-login SylabsCloud cloud.sylabs.io
  apptainer remote use SylabsCloud
  </code>

  - create container
  <code>
  export APPTAINER_CACHEDIR=/path/to/tmp
  apptainer pull --dir 'path/to/cache/dir' container_name.sif library://timniem/rna_outliers/test:sha256.7e228258c297dd9d4dd339f824f3d38845b69745285f17b5facd7052b6822781
  </code>

  - or create the container using the definitation file "rna_outlier_pipeline.def"


### Parameters

Below you'll find an explanation for the configurable parameters in nextflow.config. Before running the pipeline, make sure these parameters are set correctly:

 - singularity.cacheDir: specify path where cache files created by the pipeline will be stored
 - process.executor: use "slurm" to use slurm workload manager or "local" when you want to run the pipeline locally
 - process.container: singularity container with OUTRIDER, FRASER, MAE (including dependencies). Can be downloaded from library://timniem/rna_outliers/test:hg19hg38
 - params.genomeReferenceBuild: specify genome build used for alignment with "hg19" or "hg38"
 - params.genes_gtf: specify path to used gtf file. "gencode.v29.annotation.gtf" is used in manuscript and downloaded from GENCODE (https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_29/)
 - params.extcounts.folder: specify path to directory with external counts. A directory has to be specified. Examples of external counts can be downloaded from https://zenodo.org/records/5638707 (counts based on alignment to hg19) and https://zenodo.org/records/6078397 (counts based on alignment to hg38) (see "References"). When no external counts are used, keep "params.extcounts.amount_outrider" equal to 0 and "params.extcounts.amount_fraser" equal to 0
 - params.samplesheet: specify path to the used samplesheet (see "example_samplesheet.tsv")
 - params.fasta: specify path to used reference assembly. "GRCh37.primary_assembly.genome.fa" is used in manuscript and downloaded from GENCODE (https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_29/)
 - output: specify path to output directory. This directory specifies where the output from the RNA outlier pipeline are stored

## How to run RNA outlier pipeline

After making sure all the dependencies are installed correctly and all requirements are met, you can run the pipeline using the following instructions:

1. Set the parameters in the config in accordence to the above description
2. Run the pipeline with
```
cd path/to/nextflow
nextflow run main.nf --output "path/to/output" --samplesheet "path/to/samplesheet.tsv"
```

## How to generate interactive report
```
python embed_data.py -or sample_id_result_table_outrider.tsv -fr sample_id_result_table_fraser.tsv -t template.html -s sample_id -o sample_id.html

-or: specify outrider output .tsv of sample you want to create the report for
-fr: specify fraser output .tsv of sample you want to create the report for
-t: provide location of template html file
-s: provide sample_id
-o: provide name of report
```

## References

Brechtmann F, Mertes C, Matusevičiūtė A, et al. OUTRIDER: A Statistical Method for Detecting Aberrantly Expressed Genes in RNA Sequencing Data. Am J Hum Genet. 2018;103(6):907-917. https://doi.org/10.1016/j.ajhg.2018.10.025

Scheller, I.F., Lutz, K., Mertes, C et al. Improved detection of aberrant splicing with FRASER 2.0 and the intron Jaccard index. Am Jrnl Hum Genet 110, 12 (2023). https://doi.org/10.1016/j.ajhg.2023.10.014

Yepez, V. A., Gusic, M., Kopajtich, R., Meitinger, T., Gagneur, J., & Prokisch, H. (2021). Gene expression and splicing counts from the Yepez, Gusic et al study - non strand-specific (1.2) [Data set]. Zenodo. https://doi.org/10.5281/zenodo.7126296

Yépez, V. A., Smith, N. H., Mertes, C., & Gagneur, J. (2022). Gene expression and splicing counts from 49 tissues from GTEx v8 genome build hg38 - non-strand specific (1.0) [Data set]. Zenodo. https://doi.org/10.5281/zenodo.6078397


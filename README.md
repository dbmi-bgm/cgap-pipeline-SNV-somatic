<img src="https://github.com/dbmi-bgm/cgap-pipeline/blob/master/docs/images/cgap_logo.png" width="200" align="right">

# CGAP Pipeline for Somatic Single-Nucleotide Variants and small INDELs

This repository contains components for the CGAP pipeline for somatic single-nucleotide variants (SNVs) and small INDELs:

  * CWL workflows
  * CGAP Portal Workflows and MetaWorkflows objects
  * ECR (Docker) source files, which allow for creation of public Docker images (using `docker build`) or private dynamically-generated ECR images (using [*cgap pipeline utils*](https://github.com/dbmi-bgm/cgap-pipeline-utils/) `deploy_pipeline`)

The pipeline starts from tumor/normal `vcf` files and produces filtered and annotated `vcf` files containing SNVs and small INDELs as output. For more details check [*documentation*](https://cgap-pipeline-main.readthedocs.io/en/latest/Pipelines/Downstream/SNV_somatic/index-SNV_somatic.html "SNV somatic").

<img src="https://github.com/dbmi-bgm/cgap-pipeline/blob/master/docs/images/cgap_logo.png" width="200" align="right">

# CGAP Pipeline for Somatic Single Nucleotide Variants and small INDELs

This repository contains components for the CGAP pipeline for single nucleotide variants (SNVs) and small INDELs in somatic data:

  * CWL workflow descriptions
  * CGAP Portal *Workflow* and *MetaWorkflow* objects
  * CGAP Portal *Software*, *FileFormat*, and *FileReference* objects
  * ECR (Docker) source files, which allow for creation of public Docker images (using `docker build`) or private dynamically-generated ECR images (using [*cgap pipeline utils*](https://github.com/dbmi-bgm/cgap-pipeline-utils/) `pipeline_deploy`)

The pipeline starts from tumor/normal `vcf` files and produces filtered and annotated `vcf` files containing SNVs and small INDELs as output.
It also prioritizes and annotates putative driver mutations applying the decision tree developed by Hartwig Medical Foundation.
Documentation for all CGAP Pipelines can now be found here:
https://cgap-pipeline-main.readthedocs.io/en/latest/

#!/usr/bin/env cwl-runner

cwlVersion: v1.0

class: CommandLineTool

requirements:
  - class: InlineJavascriptRequirement

hints:
  - class: DockerRequirement
    dockerPull: ACCOUNT/driver_catalog:VERSION

baseCommand: [python3, /usr/local/bin/somatic_annot.py, driverCatalogVCF ] 

inputs:
  - id: input_vcf
    type: File
    inputBinding:
      prefix: -i
    doc: expect the path to the somatic vcf gz file

  - id: gene_panel
    type: File
    inputBinding:
      prefix: -g
    doc: gene panel configuration created by Hartwig Medical Foundation

  - id: somatic_hotspot
    type: File
    inputBinding:
        prefix: -s
    doc: vcf containing hotspot mutations

  - id: output
    type: string
    inputBinding:
        prefix: -o
    default: output.vcf
    doc: vcf file containing annotated putative drivers

outputs:

  - id: putative_drivers_vcf
    type: File
    outputBinding:
      glob: $(inputs.output).gz
    secondaryFiles:
      - .tbi

doc: |
  run driverCatalogVCF to annotate putative driver genes in input VCF file
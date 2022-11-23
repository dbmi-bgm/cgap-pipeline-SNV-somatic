#!/usr/bin/env cwl-runner

cwlVersion: v1.0

class: CommandLineTool

requirements:
  - class: InlineJavascriptRequirement

hints:
  - class: DockerRequirement
    dockerPull: ACCOUNT/hmf:VERSION

baseCommand: [python3, /usr/local/bin/somatic_annot.py, driverCatalogVCF ] 

inputs:
  - id: input
    type: File
    inputBinding:
      prefix: -i
    doc: expect the path to the somatic vcf gz file

  - id: gene_panel
    type: File
    inputBinding:
      prefix: -g
    doc: 

  - id: somatic_hotspot
    type: File
    inputBinding:
        prefix: -s

  - id: output
    type: string
    inputBinding:
        prefix: -o
    default: output.vcf

outputs:

  - id: driver_annot_vcf
    type: File
    outputBinding:
      glob: $(inputs.output).gz
    secondaryFiles:
      - .tbi
doc: |
  todo
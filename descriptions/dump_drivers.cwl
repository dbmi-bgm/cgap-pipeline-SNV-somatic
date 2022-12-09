#!/usr/bin/env cwl-runner

cwlVersion: v1.0

class: CommandLineTool

requirements:
  - class: InlineJavascriptRequirement

hints:
  - class: DockerRequirement
    dockerPull: ACCOUNT/driver_catalog:VERSION

baseCommand: [python3, /usr/local/bin/somatic_annot.py, dumpJSON ] 

inputs:
  - id: vcf
    type: File
    inputBinding:
      prefix: -v
    doc: expect the path to the somatic vcf gz file

  - id: cnv
    type: File
    inputBinding:
      prefix: -c
    doc: SV file containing somatic CNVs

  - id: output
    type: string
    inputBinding:
        prefix: -o
    default: 'output.json'
    doc: output file name

outputs:

  - id: putative_drivers_json
    type: File
    outputBinding:
      glob: $(inputs.output)
    secondaryFiles:
      - .tbi
doc: |
  JSON including reported putative drivers
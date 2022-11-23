#!/usr/bin/env cwl-runner

cwlVersion: v1.0

class: CommandLineTool

requirements:
  - class: InlineJavascriptRequirement

hints:
  - class: DockerRequirement
    dockerPull: ACCOUNT/snv_somatic:VERSION

baseCommand: [multiple_vcf-integrity-check.sh]

inputs:
  - id: input
    type: File
    inputBinding:
      position: 1
    doc: expect the path to the vcf gz file

outputs:
  - id: output
    type: File
    outputBinding:
      glob: integrity_check_$(inputs.input.basename.split(".")[0])

doc: |
  run a quick integrity check on the input vcf file

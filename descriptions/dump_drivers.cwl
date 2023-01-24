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
    doc: VCF file containing SNVs and INDELs annotated with putative driver mutations and hotspot locations

  - id: cnv
    type: File
    inputBinding:
      prefix: -c
    doc: TSV file containing CNVs identified as putative driver mutations

  - id: output
    type: string
    inputBinding:
        prefix: -o
    default: 'output.json'
    doc: output file name

outputs:

  - id: drivers_json
    type: File
    outputBinding:
      glob: $(inputs.output)
    secondaryFiles:
      - .tbi

doc: |
  run dumpJSON to create a JSON file with the reported putative driver genes
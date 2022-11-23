#!/usr/bin/env cwl-runner

cwlVersion: v1.0

class: CommandLineTool

requirements:
  - class: InlineJavascriptRequirement

hints:
  - class: DockerRequirement
    dockerPull: ACCOUNT/snv_somatic_splitter:VERSION

baseCommand: [python3, /usr/local/bin/somatic_sentieon_vcf_splitter.py, -o, snv, ind,  -o, bnd, --pass_only]

inputs:
  - id: input
    type: File
    inputBinding:
      prefix: -i
    doc: expect the path to the somatic vcf gz file

  - id: prefix
    type: string
    inputBinding:
      prefix: -p
    default: 'pre'
    doc: base name of the full PASS-filtered output vcf gz file


outputs:

  - id: indel_snv_pass
    type: File
    outputBinding:
      glob: $(inputs.prefix + "_ind_snv.vcf.gz")
    secondaryFiles:
      - .tbi
  - id: sv_bnd
    type: File
    outputBinding: 
      glob: $(inputs.prefix + "_bnd.vcf.gz")
    secondaryFiles:
      - .tbi
doc: |
  run somatic_sentieon_vcf_splitter.py to create 2 VCF files containing SNV/INDEL and SV

#!/usr/bin/env cwl-runner

cwlVersion: v1.0

class: CommandLineTool

requirements:
  - class: InlineJavascriptRequirement

hints:
  - class: DockerRequirement
    dockerPull: ACCOUNT/snv_somatic:VERSION

baseCommand: [python3, /usr/local/bin/somatic_sentieon_vcf_splitter.py]

inputs:
  - id: input
    type: File
    inputBinding:
      prefix: -i
    doc: expect the path to the somatic vcf gz file

  - id: full_pass
    type: string
    inputBinding:
      prefix: -f
    doc: base name of the full PASS-filtered output vcf gz file

  - id: snv_pass
    type: string
    inputBinding:
      prefix: -s
    doc: base name of the SNV PASS-filtered output vcf gz file

  - id: indel_pass
    type: string
    inputBinding:
      prefix: -d
    doc: base name of the INDEL PASS-filtered output vcf gz file

  - id: sv_pass
    type: string
    inputBinding:
      prefix: -v
    doc: base name of the SV PASS-filtered output vcf gz file

outputs:
  - id: full_pass
    type: File
    outputBinding:
      glob: $(inputs.full_pass + ".gz")
    secondaryFiles:
      - .tbi

  - id: snv_pass
    type: File
    outputBinding:
      glob: $(inputs.snv_pass + ".gz")
    secondaryFiles:
      - .tbi

  - id: indel_pass
    type: File
    outputBinding:
      glob: $(inputs.indel_pass + ".gz")
    secondaryFiles:
      - .tbi

  - id: sv_pass
    type: File
    outputBinding:
      glob: $(inputs.sv_pass + ".gz")
    secondaryFiles:
      - .tbi

doc: |
  run somatic_sentieon_vcf_splitter.py to create 4 filtered VCF files

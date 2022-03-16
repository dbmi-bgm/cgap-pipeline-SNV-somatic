cwlVersion: v1.0

class: Workflow

requirements:
  MultipleInputFeatureRequirement: {}

inputs:
  - id: input_vcf
    type: File
    doc: expect the path to the sample vcf gz file

  - id: full_pass
    type: string
    default: "full.vcf"
    doc: base name of the full PASS-filtered output vcf gz file

  - id: snv_pass
    type: string
    default: "snv.vcf"
    doc: base name of the SNV PASS-filtered output vcf gz file

  - id: indel_pass
    type: string
    default: "indel.vcf"
    doc: base name of the INDEL PASS-filtered output vcf gz file

  - id: sv_pass
    type: string
    default: "sv.vcf"
    doc: base name of the SV PASS-filtered output vcf gz file

outputs:

  full_pass_vcf:
    type: File
    outputSource: somatic_sentieon_vcf_splitter/full_pass

  snv_pass_vcf:
    type: File
    outputSource: somatic_sentieon_vcf_splitter/snv_pass

  indel_pass_vcf:
    type: File
    outputSource: somatic_sentieon_vcf_splitter/indel_pass

  sv_pass_vcf:
    type: File
    outputSource: somatic_sentieon_vcf_splitter/sv_pass

  vcf-check:
    type: File
    outputSource: integrity-check/output

steps:
  somatic_sentieon_vcf_splitter:
    run: somatic_sentieon_vcf_splitter.cwl
    in:
      input:
        source: input_vcf
      full_pass:
        source: full_pass
      snv_pass:
        source: snv_pass
      indel_pass:
        source: indel_pass
      sv_pass:
        source: sv_pass
    out: [full_pass, snv_pass, indel_pass, sv_pass]

  full_integrity-check:
    run: vcf-integrity-check-SNV-somatic.cwl
    in:
      input:
        source: somatic_sentieon_vcf_splitter/full_pass
    out: [output]

  snv_integrity-check:
    run: vcf-integrity-check-SNV-somatic.cwl
    in:
      input:
        source: somatic_sentieon_vcf_splitter/snv_pass
    out: [output]

  indel_integrity-check:
    run: vcf-integrity-check-SNV-somatic.cwl
    in:
      input:
        source: somatic_sentieon_vcf_splitter/indel_pass
    out: [output]

  sv_integrity-check:
    run: vcf-integrity-check-SNV-somatic.cwl
    in:
      input:
        source: somatic_sentieon_vcf_splitter/sv_pass
    out: [output]

doc: |
  run somatic_sentieon_vcf_splitter.py to create 4 filtered VCF files |
  run an integrity check on each output vcf gz

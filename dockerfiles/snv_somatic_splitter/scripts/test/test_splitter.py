#################################################################
#   Libraries
#################################################################
import sys, os
import pytest
import filecmp
from granite.lib import vcf_parser


splitter = __import__("somatic_sentieon_vcf_splitter")


def test_splitter(tmp_path):
    """
    Test for splitter
    """
    prefix = f"{tmp_path}/pre"
    out_ind_snv = f"{prefix}_ind_snv.vcf"
    out_bnd = f"{prefix}_bnd.vcf"

    # Variables and Run
    args = {"inputvcf": "test/files/in_splitter.vcf", "prefix": f"{prefix}", "output":[['snv', 'ind'], ['bnd']], "pass_only":False}

    splitter.main(args)
    result_snv_ind = filecmp.cmp(out_ind_snv, "test/files/out_splitter_snv_ind.vcf")
    result_bnd = filecmp.cmp(out_bnd, "test/files/out_splitter_bnd.vcf")

    # Test
    assert result_snv_ind == True
    assert result_bnd == True
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
    args = {"inputvcf": "test/files/in_splitter.vcf.gz", "prefix": f"{prefix}", "output":[['snv', 'ind'], ['bnd']], "pass_only":False}

    splitter.main(args)

    a = os.popen(f'bgzip -c -d {out_ind_snv}.gz')
    b = os.popen('bgzip -c -d "test/files/out_splitter_snv_ind.vcf.gz"')
    assert [row for row in a.read()] == [row for row in b.read()]
    

    c = os.popen(f'bgzip -c -d {out_bnd}.gz')
    d = os.popen(f'bgzip -c -d test/files/out_splitter_bnd.vcf.gz')
    assert [row for row in c.read()] == [row for row in d.read()]

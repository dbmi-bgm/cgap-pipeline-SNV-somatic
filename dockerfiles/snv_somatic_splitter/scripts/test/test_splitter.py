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

def test_splitter_sv(tmp_path):
    """
    Test for splitter SVs
    """
    prefix = f"{tmp_path}/sv_test"
    out_dup = f"{prefix}_dup.vcf"
    out_ins = f"{prefix}_ins.vcf"
    out_inv = f"{prefix}_inv.vcf"
    out_ins_inv_dup = f"{prefix}_dup_ins_inv.vcf"
    out_all_sv = f"{prefix}_all_sv.vcf"
    # Variables and Run
    args = {"inputvcf": "test/files/sv_test_full.vcf.gz", "prefix": f"{prefix}", "output":[['dup'], ['ins'], ['inv'], ['ins', 'inv', 'dup'], ['all_sv']], "pass_only":False}

    splitter.main(args)

    a = os.popen(f'bgzip -c -d {out_dup}.gz')
    b = os.popen('bgzip -c -d "test/files/sv_test_dup.vcf.gz"')
    assert [row for row in a.read()] == [row for row in b.read()]
    

    c = os.popen(f'bgzip -c -d {out_ins}.gz')
    d = os.popen(f'bgzip -c -d test/files/sv_test_ins.vcf.gz')
    assert [row for row in c.read()] == [row for row in d.read()]


    e = os.popen(f'bgzip -c -d {out_inv}.gz')
    f = os.popen(f'bgzip -c -d test/files/sv_test_inv.vcf.gz')
    assert [row for row in e.read()] == [row for row in f.read()]


    g = os.popen(f'bgzip -c -d {out_ins_inv_dup}.gz')
    h = os.popen(f'bgzip -c -d test/files/sv_test_ins_inv_dup.vcf.gz')
    assert [row for row in g.read()] == [row for row in h.read()]


    i = os.popen(f'bgzip -c -d {out_all_sv}.gz')
    k = os.popen(f'bgzip -c -d test/files/sv_test_full.vcf.gz')
    assert [row for row in i.read()] == [row for row in k.read()]


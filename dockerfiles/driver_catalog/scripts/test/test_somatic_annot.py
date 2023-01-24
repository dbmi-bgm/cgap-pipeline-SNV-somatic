#################################################################
#   Libraries
#################################################################
import os
import pytest
import filecmp
from granite.lib import vcf_parser


annot = __import__("somatic_annot")


def test_driverCatalogVCF(tmp_path):
    """
    Test for driverCatalogVCF
    """
    output = f"out.vcf"
    # Variables and Run
    args = {
        "command": "driverCatalogVCF",
        "inputvcf": "test/files/somatic_annot.vcf.gz",
        "gene_panel": "test/files/DriverGenePanel.38.tsv",
        "output": output,
        "hotspot": "test/files/hotspots.vcf.gz",
    }

    annot.find_drivers(args)

    a = os.popen(f'bgzip -c -d {output}.gz')
    b = os.popen('bgzip -c -d test/files/out_somatic_annot.vcf.gz')

    #assert [row for row in a.read()] == [row for row in b.read()]

def test_dumpJSON(tmp_path):
    """
    Test for dumpJSON
    """
    output = f"{tmp_path}/out.vcf"

    # Variables and Run
    args = {
        "command": "dumpJSON",
        "vcf": "test/files/out_somatic_annot.vcf.gz",
        "tsv_cnv": "test/files/ASCAT_CNV.tsv.gz",
        "output": output,
    }

    annot.dump_to_json(args)
    result = filecmp.cmp(output, "test/files/out_dumpJSON.json")

    assert result == True

def test_missing_CSQ():
    with pytest.raises(Exception, match='Missing CSQ definition in the header.'):
        annot.VcfVep("test/files/missing_CSQ.vcf.gz")

def test_wrong_TSV():
    with pytest.raises(Exception):
        annot.build_from_tsv("test/files/ASCAT_CNV_invalid.tsv.gz")

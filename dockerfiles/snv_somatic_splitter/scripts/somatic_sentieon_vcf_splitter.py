#!/usr/bin/env python3

################################################
#
#  Script to split Sentieon somatic vcf 3 ways
#   and filter for high confidence variants
#
################################################

################################################
#   Libraries
################################################
import sys, argparse, subprocess
from granite.lib import vcf_parser, shared_functions

SNV = 'snv'
INS = 'ins'
MNV = 'mnv'
DEL = 'del'
INV = 'inv'
BND = 'bnd'
CNV = 'cnv'
DUP = 'dup'
SV = 'sv'
MAV = 'mav'
ALL = "all"
INDEL = "ind"

################################################
#   Functions
################################################

def get_variant_type(vnt_obj):
    """get variant type of the variant
    :param vnt_obj: Variant object
    :type vnt_obj: vcf_parser.Vcf
    :returns: variant type
    :rtype: str
    """
    try:
        vnt_obj.get_tag_value("SVTYPE")
        return SV
    except Exception:
        return shared_functions.variant_type_ext(vnt_obj.REF, vnt_obj.ALT)


def main(args):
    # open sample VCF
    vcf = vcf_parser.Vcf(args["inputvcf"])
    vcf.header.definitions = vcf.header.definitions.replace(
        "\\t", "%09"
    )  # sentieon command with \t is an illegal character in the vcf

    full_file = snv_file = indel_file = None
    files = []
    prefix = args["prefix"] 
    
    full_file = open(f"{prefix}_full.vcf", "w")
    files.append(full_file)
    
    files_dict = {}
    if len(args["output"]) > 0:

        args["output"] = [sorted(list(set(v))) for v in args["output"]]

        for comb in args["output"]:
            name = "_".join(comb)
            files_dict[name] = open(f"{prefix}_{name}.vcf", "w")
            files.append(files_dict[name])

    for file in files:
        vcf.write_header(file)

    for vnt_obj in vcf.parse_variants():

        if vnt_obj.FILTER != "PASS" and args["pass_only"] == True:
            continue

        variant_type = get_variant_type(vnt_obj)

        if variant_type == SNV:
            for comb in args["output"]:
                if SNV in comb:
                    vcf.write_variant(files_dict["_".join(comb)], vnt_obj)

        if variant_type in [INS, DEL, MNV]:
            for comb in args["output"]:
                if INDEL in comb:
                    vcf.write_variant(files_dict["_".join(comb)], vnt_obj)

        #print(variant_type)
        if variant_type == SV and len(list(files_dict.keys())) > 0:
            sv_type = shared_functions.variant_type_sv(vnt_obj)
            for sv_comb in args["output"]:
                #print(sv_type, " ", sv_comb)
                if sv_type in sv_comb:
                    vcf.write_variant(files_dict["_".join(sv_comb)], vnt_obj)
                if ALL in sv_comb:
                    vcf.write_variant(files_dict["_".join(sv_comb)], vnt_obj)

    for file in files:
        file.close()
        subprocess.run(["bgzip", file.name])
        subprocess.run(["tabix", file.name + ".gz"])


################################################
#   MAIN
################################################
if __name__ == "__main__":

    parser = argparse.ArgumentParser(description="")

    parser.add_argument("-i", "--inputvcf", help="input sample vcf", required=True)
    parser.add_argument("-p", "--prefix", help="prefix for the output files", required=True)
    parser.add_argument(
        "-o",
        "--output",
        help="variants type",
        nargs="+",
        action = 'append',
        required=False,
        choices=[CNV, DEL, DUP, BND, INV, INS, DEL, ALL, SNV, INDEL],
    )

    parser.add_argument(
        "-f",
        "--pass_only",
        help="only PASS variants",
        action = 'store_true',
        required=False,
    )
    
    args = vars(parser.parse_args())

    main(args)

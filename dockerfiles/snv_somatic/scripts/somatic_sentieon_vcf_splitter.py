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
from granite.lib import vcf_parser

################################################
#   Functions
################################################

def variant_type(vnt_obj):
    ''' Returns SV, SNV, or INDEL depending on variant type '''
    try:
        vnt_obj.get_tag_value("SVTYPE")
        return("SV")
    except:
        if len(vnt_obj.REF) == 1:
            if len(vnt_obj.ALT) == 1:
                return("SNV")
            else:
                return("INDEL")
        elif len(vnt_obj.ALT) == 1:
            return("INDEL")

def variant_write(TYPE, VCF_obj, VNT_obj, full_FILE, snv_FILE, indel_FILE, sv_FILE):
    ''' Write to correct files, given variant type '''
    VCF_obj.write_variant(full_FILE, VNT_obj)

    if TYPE == "SNV":
        VCF_obj.write_variant(snv_FILE, VNT_obj)
    elif TYPE == "INDEL":
        VCF_obj.write_variant(indel_FILE, VNT_obj)
    elif TYPE == "SV":
        VCF_obj.write_variant(sv_FILE, VNT_obj)
    else:
        raise Exception('No variant type identified for '+VNT_obj.CHROM+':'+str(VNT_obj.POS)+' '+VNT_obj.REF+' '+VNT_obj.ALT)

def main(args):

    #open sample VCF
    vcf = vcf_parser.Vcf(args['inputvcf'])

    with open(args['full'], 'w') as full, open(args['snv'], 'w') as SNV, open(args['indel'], 'w') as INDEL, open(args['sv'], 'w') as SV:
        for file in [full, SNV, INDEL, SV]:
            vcf.write_header(file)

        for vnt_obj in vcf.parse_variants():
            if vnt_obj.FILTER == "PASS":
                variant_write(variant_type(vnt_obj), vcf, vnt_obj, full, SNV, INDEL, SV)

    for file in [args['full'], args['snv'], args['indel'], args['sv']]:
        subprocess.run(["bgzip", file])
        subprocess.run(["tabix", file+".gz"])

################################################
#   MAIN
################################################
if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='')

    parser.add_argument('-i', '--inputvcf', help='input sample vcf', required=True)
    parser.add_argument('-f', '--full', help='output VCF file for all variant types', required=True)
    parser.add_argument('-s', '--snv', help='output VCF file for SNVs', required=True)
    parser.add_argument('-d', '--indel', help='output VCF file for INDELs', required=True)
    parser.add_argument('-v', '--sv', help='output VCF file for SVs', required=True)

    args = vars(parser.parse_args())

    main(args)

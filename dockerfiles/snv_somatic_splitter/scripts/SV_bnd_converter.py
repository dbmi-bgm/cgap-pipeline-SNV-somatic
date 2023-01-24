#!/usr/bin/env python3

################################################
#
#  Script to convert Sentieon SVs to BEDPE
#
################################################



################################################
#   Libraries
################################################
import sys, os, argparse, subprocess
from granite.lib import vcf_parser
import re

################################################
#   Classes
################################################


class bnd:
    def __init__(self, vnt_obj):
        self.CHROM = vnt_obj.CHROM
        self.POS = vnt_obj.POS
        self.ID = vnt_obj.ID
        self.REF = vnt_obj.REF
        self.ALT = vnt_obj.ALT.split(",")
        self.QUAL = vnt_obj.QUAL
        self.INFO = vnt_obj.INFO
        self.STRAND = [ "-" if any(alt.startswith(x) for x in ["[", "]"]) else "+" for alt in vnt_obj.ALT.split(",") ]
        self.MATES = mate_lister(vnt_obj)


#vcf = vcf_parser.Vcf("/Users/phil_hms/Desktop/somatic/SV_test_tumor_normal_with_panel.vcf")
# vcf = vcf_parser.Vcf("/Users/phil_hms/Desktop/somatic/multi_match.vcf")
# for vnt_obj in vcf.parse_variants():
#     if vnt_obj.get_tag_value("SVTYPE") == "BND":
#         bnd_obj = bnd(vnt_obj)
#         print(bnd_obj.ID,bnd_obj.STRAND,bnd_obj.MATES)


################################################
#   Functions
################################################


def mate_lister(vnt_obj):
    try:
        return vnt_obj.get_tag_value("MATEID").split(",")
    except:
        return []

def mate_matcher(bnd_obj, mate_dict):
    if bnd_obj.ID not in mate_dict:
        if len(bnd_obj.MATES) == 1:
            mate_dict[bnd_obj.ID] = {bnd_obj.MATES[0]:bnd_obj}

        else:
            raise ValueError('\nERROR: no MATEID for '+bnd_obj.ID+'. Exiting\n')
    else:
        raise ValueError('\nERROR: duplicated bnd for '+bnd_obj.ID+'. Exiting\n')

def vcf_scanner(vcf_input):
    #vcf = vcf_parser.Vcf(args['inputvcf'])
    vcf = vcf_parser.Vcf(vcf_input)
    mate_dict = {}
    bnd_list = []

    for vnt_obj in vcf.parse_variants():
        if vnt_obj.get_tag_value("SVTYPE") == "BND":
            bnd_list.append(vnt_obj.ID)
            bnd_obj = bnd(vnt_obj)
            mate_matcher(bnd_obj,mate_dict)
        #more robust?
        else:
            #so far, I have only seen SVTYPE=INS in tumor only analyses
            #we will only support paired breakpoints at this point with BEDPE
            if vnt_obj.get_tag_value("SVTYPE") == "INS":
                pass

    return mate_dict, bnd_list

#create list of IDs and run the list removing keys from the dict.
#dictionary.get? to check if key is there...
def mate_checker(mateID, mate_dict):
    # mate_dict has {bnd1: {bnd2: <bnd class object for bnd1>}}
    # want to check that bnd1's mate (bnd2) has bnd1 as its mate in its entry in the dictionary
    # only want to do this once for each pair, so finished_list is imporant there
    mate_match = list(mate_dict[mateID].keys())[0]
    if mateID == list(mate_dict[mate_match].keys())[0]:
        # if they do match reciprocally, we want to store each of their vnt_objs from the original vcf file
        pair1 = list(mate_dict[mateID].values())[0]
        pair2 = list(mate_dict[mate_match].values())[0]

        # add the second mate to the finished_dict so that we only do the work once
        #finished_list.append(mate_match)

        #remove the second mate from the mate_dict to prevent doubling the work
        del mate_dict[mate_match]

        #then return the pair of vnt_objs
        return pair1, pair2, mate_dict

    else:
        raise ValueError('\nERROR: bnd mates not reciprocal for '+mateID+'. Exiting\n')

def strand_finder(bnd):
    #if alt.startswith [ or ]
    alt = re.findall(r'([])([])', bnd.ALT)
    #print(alt)
    # bnd square brackets must be in the same orientation
    try:
        assert alt[0] == alt[1], 'ALT square brackets must be in same orientation: [[ or ]]'

        # + strand if breakpoint after position (e.g., A[), and - strand if breakpoint before position (e.g, [)
        if bnd.ALT.startswith(alt[0]):
            return('-')
        else:
            return('+')

    except AssertionError as e:
        raise

def create_bedpe(pair1, pair2):
    chromosomes = ["chr1",
                    "chr2",
                    "chr3",
                    "chr4",
                    "chr5",
                    "chr6",
                    "chr7",
                    "chr8",
                    "chr9",
                    "chr10",
                    "chr11",
                    "chr12",
                    "chr13",
                    "chr14",
                    "chr15",
                    "chr16",
                    "chr17",
                    "chr18",
                    "chr19",
                    "chr20",
                    "chr21",
                    "chr22",
                    "chrX",
                    "chrY"
                ]

    #only want the main chromosomes
    if pair1.CHROM not in chromosomes or pair2.CHROM not in chromosomes:
        pass
    else:
        CHROM1 = pair1.CHROM
        START1 = str(int(pair1.POS)-1)
        END1 = str(pair1.POS)
        CHROM2 = pair2.CHROM
        START2 = str(int(pair2.POS)-1)
        END2 = str(pair2.POS)
        NAME = pair1.ID
        SCORE = str(pair1.QUAL)
        STRAND1 = pair1.STRAND[0]
        STRAND2 = pair2.STRAND[0]

        if CHROM1 == CHROM2:
            if STRAND1 == "+":
                if STRAND2 == "-":
                    SVTYPE = "DEL" # DEL = +-
                else:
                    SVTYPE = "h2hINV" # h2hINV = ++

            else:
                if STRAND2 == "+":
                    SVTYPE = "DUP" # DUP = -+
                else:
                    SVTYPE = "t2tINV" # t2tINV = --

        else:
            SVTYPE = "TRA"

        bedpe_variant = [CHROM1, START1, END1, CHROM2, START2, END2, NAME, SCORE, STRAND1, STRAND2, SVTYPE, 'Sentieon_TNscope']
        return bedpe_variant

def main(args):

    # fill the mate_dict with pairs
    mate_dict, bnd_list = vcf_scanner(args['inputvcf'])

    # get matched pairs from mate_dict
    #global finished_list
    #finished_list = []

    #with open(out_file_name, 'w') as fo:
    with open(args['outputBEDPE'], 'w') as fo:
        fo.write('\t'.join(["chrom1", "start1","end1","chrom2","start2","end2","sv_id","pe_support","strand1","strand2","svclass","svmethod"])+'\n')
        for variant in bnd_list:
            pair1 = pair2 = False
            if variant in mate_dict:
                #as we move down the list we will have variants that have been moved to finished_list
                pair1, pair2, mate_dict = mate_checker(variant, mate_dict)
                #print(pair1.CHROM, pair2.CHROM)


            if pair1 and pair2:
                bedpe_variant = create_bedpe(pair1, pair2)

                if bedpe_variant:
                    fo.write('\t'.join(bedpe_variant)+'\n')

################################################
#   MAIN
################################################
if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='')

    parser.add_argument('-i', '--inputvcf',  help='input SV vcf', required=True)
    parser.add_argument('-o', '--outputBEDPE', help='output bedPE', required=True)

    args = vars(parser.parse_args())

    main(args)

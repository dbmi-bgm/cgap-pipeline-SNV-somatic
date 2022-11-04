
################################################
#
#  Script to TODO
#
################################################



################################################
#   Libraries
################################################
import sys, os, argparse, subprocess
from granite.lib import vcf_parser
import re
import pandas


def create_hotspot_key(vnt_obj):
    return f'{vnt_obj.CHROM}_{vnt_obj.POS}_{vnt_obj.REF}_{vnt_obj.ALT}'


def main(args):

    vcf_obj = vcf_parser.Vcf(args['inputvcf'])
    hotspot_vcf = vcf_parser.Vcf(args['hotspot'])
    
    hotspot_dict = {}
    results = []
    vcf_obj.header.add_tag_definition('##INFO=<ID=DRIVER,Number=1,Type=Integer,Description="Flag if variant is a putative driver mutation calculated based on Hartwig\'s decision tree (1 if driver, otherwise 0)">', tag_type="INFO")

    for vnt_obj in hotspot_vcf.parse_variants():
        hotspot_dict[create_hotspot_key(vnt_obj)] = vnt_obj
    
    

    CSQ_tag  = 'CSQ'
    gene_field= 'SYMBOL'
    consequence_field = 'Consequence'


    idx_gene = vcf_obj.header.get_tag_field_idx(CSQ_tag, gene_field) 
    idx_variant_cons = vcf_obj.header.get_tag_field_idx(CSQ_tag, consequence_field) 


    gene_panel = pandas.read_csv(args['gene_panel'], sep ='\t')
    report_missense  = gene_panel.loc[gene_panel["reportMissense"] == True , ["gene"]]
    report_nonsense  = gene_panel.loc[gene_panel["reportNonsense"] == True , ["gene"]]
    with open(args["output"], "w") as output:

        vcf_obj.write_header(output)
        for vnt_obj in vcf_obj.parse_variants():

            genes = []
            transcripts = vnt_obj.get_tag_value(CSQ_tag).split(",")
            gene_symbol_missense = []
            gene_symbol_nonsense = []
            hotspot = False

            hotspot_key = create_hotspot_key(vnt_obj)
            hotspot_candidate = True if hotspot_key in hotspot_dict.keys() else False
            gene_hostpot = None

            if hotspot_candidate:
                gene_hostpot = hotspot_dict[hotspot_key].get_tag_value("input").split("|")[0]

            DRIVER = 0 
            for transcript in transcripts:
                
                transcript_fields = transcript.split("|")
                gene = transcript_fields[idx_gene]
                consequence = transcript_fields[idx_variant_cons]
                

                if ('missense_variant' in consequence and ((report_missense['gene'] == gene)).any()) or (('frameshift_variant' in consequence or 'stop_gained' in consequence) and ((report_nonsense['gene'] == gene)).any()):
                    DRIVER = 1
                    break


            vnt_obj.add_tag_info(f"DRIVER={DRIVER}")
            vcf_obj.write_variant(output, vnt_obj)


################################################
#   MAIN
################################################
if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="")

    parser.add_argument("-i", "--inputvcf", help="input vcf", required=True)
    parser.add_argument("-g", "--gene_panel", help="gene panel", required=True)
    parser.add_argument(
        "-o",
        "--output",
        help="output VCF with annotated putative drivers",
        required=True,
    )

    parser.add_argument(
        "-c",
        "--cgap_genes",
        help="table of genes supported by CGAP",
        required=True,
    )

    parser.add_argument(
        "-s",
        "--hotspot",
        help="",
        required=True,
    )
    
    args = vars(parser.parse_args())

    main(args)



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
import json

def create_hotspot_key(vnt_obj):
    return f'{vnt_obj.CHROM}_{vnt_obj.POS}_{vnt_obj.REF}_{vnt_obj.ALT}'

def create_record(vcf_obj, gene):
    return {
        "chrom": vcf_obj.CHROM,
        "pos": vcf_obj.POS,
        "ref": vcf_obj.REF,
        "alt": vcf_obj.ALT,
        "gene": gene}
        
def main(args):

    vcf_obj = vcf_parser.Vcf(args['inputvcf'])
    hotspot_vcf = vcf_parser.Vcf(args['hotspot'])
    
    hotspot_dict = {}
    results = []
    for vnt_obj in hotspot_vcf.parse_variants():
        hotspot_dict[create_hotspot_key(vnt_obj)] = vnt_obj
    
    

    CSQ_tag  = 'CSQ'
    gene_field= 'SYMBOL'
    consequence_field = 'Consequence'
    cDNA_field = 'HGVSc'
    protein_position_field = 'HGVSp'

    idx_gene = vcf_obj.header.get_tag_field_idx(CSQ_tag, gene_field) 
    idx_variant_cons = vcf_obj.header.get_tag_field_idx(CSQ_tag, consequence_field) 
    idx_cDNA_field = vcf_obj.header.get_tag_field_idx(CSQ_tag, cDNA_field)
    idx_protein_position_field = vcf_obj.header.get_tag_field_idx(CSQ_tag, protein_position_field) 

    gene_panel = pandas.read_csv(args['gene_panel'], sep ='\t')
    report_missense  = gene_panel.loc[gene_panel["reportMissense"] == True , ["gene"]].tolist()
    report_nonsense  = gene_panel.loc[gene_panel["reportNonsense"] == True , ["gene"]].tolist()

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


        for transcript in transcripts:
            transcript_fields = transcript.split("|")
            gene = transcript_fields[idx_gene]
            consequence = transcript_fields[idx_variant_cons]
            cDNA = transcript_fields[idx_cDNA_field]
            protein_position = transcript_fields[idx_protein_position_field]
            

            if 'missense_variant' in consequence and ((report_missense['gene'] == gene)).any():
                gene_symbol = report_missense[(report_missense['gene'] == gene )]["gene"].tolist()[0]
                if gene_symbol not in genes and gene_hostpot != gene_symbol :
                    genes += [gene_symbol]
    
        
            if ('frameshift_variant' in consequence or 'stop_gained' in consequence) and ((report_nonsense['gene'] == gene)).any():
                gene_symbol = report_nonsense[(report_nonsense['gene'] == gene )]["gene"].tolist()[0]
                if gene_symbol not in genes and gene_hostpot != gene_symbol:
                    genes += [gene_symbol]
                
        
        
        for gene in genes:

            results.append(create_record(vnt_obj, gene))



    json_object = json.dumps(results, indent=4)
    
    with open("sample.json", "w") as outfile:
        outfile.write(json_object)            


            

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
        help="driver catalog",
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


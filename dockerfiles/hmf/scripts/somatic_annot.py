import csv
import json
import argparse
from granite.lib.vcf_parser import Vcf
import pandas
import subprocess

gene_field= 'SYMBOL'
CSQ_tag  = 'CSQ'
consequence_field = 'Consequence'
HGVSc_field = 'HGVSc' 
HGVSp_field = 'HGVSp'
CANONICAL_field = "CANONICAL"
DRIVER_field = "DRIVER"
DRIVER_DEFAULT = "-"
CATEGORY_SNV_INDEL = "mutational"

class Driver:
    def __init__(self, gene, category=None):
        self.category = category
        self.gene = gene

    def to_dict(self):
        pass

class DriverSnvIndel(Driver):
    def __init__(self,vnt_obj, HGVSc=None, HGVSp=None, allele_fraction=None, mutation_type=None, **kwargs):

        self.chrom = vnt_obj.CHROM
        self.pos = vnt_obj.POS
        self.ref = vnt_obj.REF
        self.alt = vnt_obj.ALT
        self.transcript_consequence = HGVSc.split(":")[1]
        self.protein_mutation = HGVSp.split(":")[1]
        self.ens_gene = HGVSc.split(":")[0]
        self.ens_st = HGVSp.split(":")[0]
        self.allele_fraction = allele_fraction
        self.mutation_type = mutation_type 
        super().__init__(**kwargs)

    def to_dict(self):
        fields = vars(self)
        return {k: v for (k,v) in fields.items() if k not in ['ens_gene', 'ens_st']}

    def __str__(self):
        return f"{self.gene}|{self.ens_gene}|{self.ens_st}"

class DriverCnv(Driver):
    def __init__(self, chrom, start, end, **kwargs):
        self.chr = chrom
        self.start = start
        self.end = end
        super().__init__(**kwargs)

    def to_dict(self):
        return vars(self)



class HartwigDecisionTree:
    def __init__(self, driver_panel, hotspot_mutations):
        gene_panel = pandas.read_csv(driver_panel, sep ='\t')
        self.report_missense = [gene[0] for gene in gene_panel.loc[gene_panel["reportMissense"] == True , ["gene"]].values.tolist()]
        self.report_nonsense = [gene[0] for gene in  gene_panel.loc[gene_panel["reportNonsense"] == True , ["gene"]].values.tolist()]
        self.hotspot_dict = {}
        hotspot_vcf = Vcf(hotspot_mutations)

        for vnt_obj in hotspot_vcf.parse_variants():
            self.hotspot_dict[self.__create_hotspot_key(vnt_obj)] = vnt_obj
        

    def __create_record_snv_indel(self, vcf_obj, record):
        variant = {
            "chrom": vcf_obj.CHROM,
            "pos": vcf_obj.POS,
            "ref": vcf_obj.REF,
            "alt": vcf_obj.ALT}
        variant.update(record)
        return variant

    def __create_hotspot_key(self, vnt_obj):
        return f'{vnt_obj.CHROM}_{vnt_obj.POS}_{vnt_obj.REF}_{vnt_obj.ALT}'

    def build(self,input_vcf, save_vcf = None):
        vcf_obj = Vcf(input_vcf)
        all_drivers = []


        idx_gene = vcf_obj.header.get_tag_field_idx(CSQ_tag, gene_field)
        idx_variant_cons = vcf_obj.header.get_tag_field_idx(CSQ_tag, consequence_field)
        idx_HGVSc =  vcf_obj.header.get_tag_field_idx(CSQ_tag, HGVSc_field)
        idx_HGVSp = vcf_obj.header.get_tag_field_idx(CSQ_tag, HGVSp_field)
        idx_CANONICAL = vcf_obj.header.get_tag_field_idx(CSQ_tag, CANONICAL_field)

        def query( vnt_obj):
            
            transcripts = vnt_obj.get_tag_value(CSQ_tag).split(",")
            drivers = []
            for transcript in transcripts:
                transcript_fields = transcript.split("|")
                gene = transcript_fields[idx_gene]
                consequence = transcript_fields[idx_variant_cons]
                HGVSc = transcript_fields[idx_HGVSc]
                HGVSp = transcript_fields[idx_HGVSp]
                CANONICAL = transcript_fields[idx_CANONICAL]

                if ((('missense_variant' in consequence) and (gene in self.report_missense)) or (('frameshift_variant' in consequence or 'stop_gained' in consequence) and  (gene in self.report_nonsense ))) and CANONICAL == "YES":
                    
                        tumor_sample = vnt_obj.IDs_genotypes[0]
                        af = vnt_obj.get_genotype_value(tumor_sample, 'AF')
                        drivers.append(DriverSnvIndel(vnt_obj = vnt_obj, HGVSc = HGVSc, HGVSp = HGVSp, allele_fraction = af,  mutation_type = consequence, category = CATEGORY_SNV_INDEL, gene = gene))
            return drivers

        if save_vcf == None:
            for vnt_obj in vcf_obj.parse_variants():
                all_drivers += query(vnt_obj)
        else:
            with open(save_vcf, "w") as output:
                vcf_obj.header.add_tag_definition(f'##INFO=<ID={DRIVER_field},Number=1,Type=Integer,Description="Flag if variant is a putative driver mutation calculated based on Hartwig\'s decision tree (1 if driver, otherwise 0)">', tag_type="INFO")
                vcf_obj.write_header(output)
                for vnt_obj in vcf_obj.parse_variants():
                    DRIVER = DRIVER_DEFAULT
                    drivers = query(vnt_obj)
                    all_drivers += drivers
                    if len(drivers) > 0:
                        DRIVER = ",".join([str(d) for d in drivers])
                    vnt_obj.add_tag_info(f"{DRIVER_field}={DRIVER}")
                    vcf_obj.write_variant(output, vnt_obj)
            

            subprocess.run(["bgzip", save_vcf])
            subprocess.run(["tabix", save_vcf+".gz"])


        return all_drivers


def build_from_tsv(input_tsv):
    input_file = csv.DictReader(open(input_tsv),  delimiter = '\t')
    results = []
    fields = ["chrom", "start", "end", "gene", "category"]
    for inn in input_file:
        record = dict(inn)
        if sorted(record.keys()) != sorted(fields):
            raise Exception(f"Wrong fields in the CNV file expected {fields}, got {list(record.keys())}")
        results += [DriverCnv(**record).to_dict()]
    return results

def build_from_vcf(input_vcf):
    records  = []
    vcf_obj = Vcf(input_vcf)
    all_drivers = []
    idx_gene = vcf_obj.header.get_tag_field_idx(CSQ_tag, gene_field)
    idx_variant_cons = vcf_obj.header.get_tag_field_idx(CSQ_tag, consequence_field)
    idx_HGVSc =  vcf_obj.header.get_tag_field_idx(CSQ_tag, HGVSc_field)
    idx_HGVSp = vcf_obj.header.get_tag_field_idx(CSQ_tag, HGVSp_field)
    idx_CANONICAL = vcf_obj.header.get_tag_field_idx(CSQ_tag, CANONICAL_field)

    for vnt_obj in vcf_obj.parse_variants():
        driver = vnt_obj.get_tag_value(DRIVER_field)

        if driver != DRIVER_DEFAULT:
            drivers = driver.split(",")
            transcripts = vnt_obj.get_tag_value(CSQ_tag).split(",")

            for transcript in transcripts:
                transcript_fields = transcript.split("|")
                gene = transcript_fields[idx_gene]
                ens_gene = transcript_fields[idx_HGVSc].split(":")[0]
                ens_st = transcript_fields[idx_HGVSp].split(":")[0]

                candidate = f"{gene}|{ens_gene}|{ens_st}"

                if candidate in drivers:
                    consequence = transcript_fields[idx_variant_cons]
                    HGVSc = transcript_fields[idx_HGVSc]
                    HGVSp = transcript_fields[idx_HGVSp]
                    tumor_sample = vnt_obj.IDs_genotypes[0]
                    af = vnt_obj.get_genotype_value(tumor_sample, 'AF')
                    records.append(DriverSnvIndel(vnt_obj = vnt_obj, HGVSc = HGVSc, HGVSp = HGVSp, allele_fraction = af,  mutation_type = consequence, category = CATEGORY_SNV_INDEL, gene = gene).to_dict())


    return records
                    


def find_drivers(args):
    hartwig = HartwigDecisionTree(args['gene_panel'], args["hotspot"])
    hartwig.build(args["inputvcf"], args["output"])



def dump_to_json(args):
    results = build_from_vcf(args["vcf"])
    cnvs  = build_from_tsv(args["tsv_cnv"])
    results += cnvs
    json_object = json.dumps(results, indent=4)
    # Writing to sample.json
    with open(args["output"], "w") as outfile:
        outfile.write(json_object)

    


################################################
#   MAIN
################################################
if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="")
    subparsers = parser.add_subparsers(help='sub-command help',  dest='command')
    driver_catalog = subparsers.add_parser('driverCatalogVCF', help='a help')

    driver_catalog.add_argument("-i", "--inputvcf", help="input vcf", required=True)
    driver_catalog.add_argument("-g", "--gene_panel", help="gene panel", required=True)
    driver_catalog.add_argument(
        "-o",
        "--output",
        help="output VCF with annotated putative drivers",
        required=True,
    )


    driver_catalog.add_argument(
        "-s",
        "--hotspot",
        help="VCF containing hotspot mutations",
        required=True,
    )


    dump_json = subparsers.add_parser('dumpJSON', help='a help')

    dump_json.add_argument("-v", "--vcf", help="input vcf", required=True)
    dump_json.add_argument("-c", "--tsv_cnv", help="input vcf", required=True)
    dump_json.add_argument(
        "-o",
        "--output",
        help="output VCF with annotated putative drivers",
        required=True,
    )

    args = parser.parse_args()
    if args.command == "driverCatalogVCF":
        find_drivers(vars(args))
    
    if args.command == "dumpJSON":
        dump_to_json(vars(args))



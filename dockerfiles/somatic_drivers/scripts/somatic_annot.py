import csv
import json
import argparse
from granite.lib.vcf_parser import Vcf
import pandas
import subprocess
import gzip

CSQ_tag = "CSQ"
gene_field = "SYMBOL"
CONSEQUENCE_field = "Consequence"
HGVSc_field = "HGVSc"
HGVSp_field = "HGVSp"
CANONICAL_field = "CANONICAL"
DRIVER_field = "DRIVER"
HOTSPOT_field = "HOTSPOT"
CATEGORY_SNV_INDEL = "mutational"



class Driver:
    def __init__(self, gene, category=None):
        self.category = category
        self.gene = gene

    def to_dict(self):
        pass


class DriverSnvIndel(Driver):
    def __init__(
        self,
        vnt_obj,
        HGVSc=None,
        HGVSp=None,
        allele_fraction=None,
        mutation_type=None,
        hotspot = False,
        **kwargs
    ):

        self.chrom = vnt_obj.CHROM
        self.pos = vnt_obj.POS
        self.ref = vnt_obj.REF
        self.alt = vnt_obj.ALT
        self.transcript_consequence = HGVSc.split(":")[1] if HGVSc else None
        self.protein_mutation = HGVSp.split(":")[1] if HGVSp else None
        self.ens_gene = HGVSc.split(":")[0] if HGVSc else None
        self.ens_st = HGVSp.split(":")[0] if HGVSp else None
        self.allele_fraction = allele_fraction
        self.mutation_type = mutation_type
        self.hotspot = hotspot
        super().__init__(**kwargs)

    def to_dict(self):
        fields = vars(self)
        return {k: v for (k, v) in fields.items() if k not in ["ens_gene", "ens_st", "hotspot"]}

    def __str__(self):
        return f"{self.gene}|{self.ens_gene}|{self.ens_st}"


class DriverCnv(Driver):
    def __init__(self, chr, start, end, **kwargs):
        self.chr = chr
        self.start = start
        self.end = end
        super().__init__(**kwargs)

    def to_dict(self):
        return vars(self)


class HartwigDecisionTree:
    def __init__(self, driver_panel, hotspot_mutations):
        gene_panel = pandas.read_csv(driver_panel, sep="\t")
        self.report_missense = [
            gene[0]
            for gene in gene_panel.loc[
                gene_panel["reportMissense"] == True, ["gene"]
            ].values.tolist()
        ]
        self.report_nonsense = [
            gene[0]
            for gene in gene_panel.loc[
                gene_panel["reportNonsense"] == True, ["gene"]
            ].values.tolist()
        ]

        self.non_canonical = gene_panel.set_index('gene')['additionalReportedTranscripts'].dropna().to_dict()
        self.hotspot_dict = {}
        hotspot_vcf = Vcf(hotspot_mutations)

        for vnt_obj in hotspot_vcf.parse_variants():
            self.hotspot_dict[self.__create_hotspot_key(vnt_obj)] = vnt_obj

    def __create_record_snv_indel(self, vcf_obj, record):
        variant = {
            "chrom": vcf_obj.CHROM,
            "pos": vcf_obj.POS,
            "ref": vcf_obj.REF,
            "alt": vcf_obj.ALT,
        }
        variant.update(record)
        return variant

    def __create_hotspot_key(self, vnt_obj):
        return f"{vnt_obj.CHROM}_{vnt_obj.POS}_{vnt_obj.REF}_{vnt_obj.ALT}"

    def build(self, input_vcf, save_vcf=None):
        vcf_obj = Vcf(input_vcf)
        all_drivers = []

        idx_gene = vcf_obj.header.get_tag_field_idx(CSQ_tag, gene_field)
        idx_variant_cons = vcf_obj.header.get_tag_field_idx(CSQ_tag, CONSEQUENCE_field)
        idx_HGVSc = vcf_obj.header.get_tag_field_idx(CSQ_tag, HGVSc_field)
        idx_HGVSp = vcf_obj.header.get_tag_field_idx(CSQ_tag, HGVSp_field)
        idx_CANONICAL = vcf_obj.header.get_tag_field_idx(CSQ_tag, CANONICAL_field)

        def check_non_canonical(gene, transcript):
            if gene in self.non_canonical.keys():
                if transcript in self.non_canonical[gene]:
                    return True
            return False
        
        def check_hotspot(vnt_obj):
            vnt_key = self.__create_hotspot_key(vnt_obj)
            if vnt_key in self.hotspot_dict.keys():
                return self.hotspot_dict[vnt_key]
            else:
                return None      
        def query_hotspot(vnt_obj):
            hotspot = check_hotspot(vnt_obj)
            if hotspot:
                return hotspot.get_tag_value("input")
            else: 
                return None
            
        def query(vnt_obj):
            transcripts = vnt_obj.get_tag_value(CSQ_tag).split(",")
            drivers = []

            

            for transcript in transcripts:
                transcript_fields = transcript.split("|")
                gene = transcript_fields[idx_gene]
                #HGVSc = transcript_fields[idx_HGVSc]
                #transcript = HGVSc.split(":")[0]                    
                consequence = transcript_fields[idx_variant_cons]
                #HGVSp = transcript_fields[idx_HGVSp]
                ADD_TRANSCRIPT = False
                if transcript_fields[idx_CANONICAL] == "YES":
                    ADD_TRANSCRIPT = True
                else:
                   ADD_TRANSCRIPT = check_non_canonical(gene, transcript)
                if (
                    (
                        ("missense_variant" in consequence)
                        and (gene in self.report_missense)
                    )
                    or (
                        (
                            "frameshift_variant" in consequence
                            or "stop_gained" in consequence
                        )
                        and (gene in self.report_nonsense)
                    )
                ) and ADD_TRANSCRIPT == True:
                    #tumor_sample = vnt_obj.IDs_genotypes[0]
                    #af = vnt_obj.get_genotype_value(tumor_sample, "AF")
                    drivers.append(
                        DriverSnvIndel(
                            vnt_obj=vnt_obj,
                            gene=gene)
                    )
            return drivers

        if save_vcf == None:
            for vnt_obj in vcf_obj.parse_variants():
                all_drivers += query(vnt_obj)

        else:
            with open(save_vcf, "w") as output:
                vcf_obj.header.add_tag_definition(
                    f'##INFO=<ID={DRIVER_field},Number=1,Type=String,Description="Putative driver mutation calculated based on Hartwig\'s decision tree. Format: ">',
                    tag_type="INFO",
                )
                vcf_obj.header.add_tag_definition(
                    f'##INFO=<ID={HOTSPOT_field},Number=1,Type=String,Description="Somatic hotspot location. Format: ">',
                    tag_type="INFO",
                )
                vcf_obj.write_header(output)
                for vnt_obj in vcf_obj.parse_variants():
                    drivers = query(vnt_obj)
                    
                    all_drivers += drivers
                    if len(drivers) > 0:
                        DRIVER = ",".join([str(d) for d in drivers])
                        vnt_obj.add_tag_info(f"{DRIVER_field}={DRIVER}")
                    
                    hotspot = query_hotspot(vnt_obj)
                    if hotspot:
                        vnt_obj.add_tag_info(f"{HOTSPOT_field}={hotspot}")
                    vcf_obj.write_variant(output, vnt_obj)

            subprocess.run(["bgzip", save_vcf])
            subprocess.run(["tabix", save_vcf + ".gz"])

        return all_drivers


def build_from_tsv(input_tsv):
    input_file = csv.DictReader(gzip.open(input_tsv, mode="rt"), delimiter="\t")
    results = []
    fields = ["chr", "start", "end", "gene", "category"]
    for inn in input_file:
        record = dict(inn)
        if sorted(record.keys()) != sorted(fields):
            raise Exception(
                f"Wrong fields in the CNV file expected {fields}, got {list(record.keys())}"
            )
        results += [DriverCnv(**record).to_dict()]
    return results


def build_from_vcf(input_vcf):
    records = []
    vcf_obj = Vcf(input_vcf)
    all_drivers = []
    idx_gene = vcf_obj.header.get_tag_field_idx(CSQ_tag, gene_field)
    idx_variant_cons = vcf_obj.header.get_tag_field_idx(CSQ_tag, CONSEQUENCE_field)
    idx_HGVSc = vcf_obj.header.get_tag_field_idx(CSQ_tag, HGVSc_field)
    idx_HGVSp = vcf_obj.header.get_tag_field_idx(CSQ_tag, HGVSp_field)
    idx_CANONICAL = vcf_obj.header.get_tag_field_idx(CSQ_tag, CANONICAL_field)

    for vnt_obj in vcf_obj.parse_variants():
        try:
            hotspot = vnt_obj.get_tag_value(HOTSPOT_field)
        try:
            driver = vnt_obj.get_tag_value(DRIVER_field)
        except Exception:
            continue 
        tumor_sample = vnt_obj.IDs_genotypes[0]
        af = vnt_obj.get_genotype_value(tumor_sample, "AF")
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

                    records.append(
                        DriverSnvIndel(
                            vnt_obj=vnt_obj,
                            HGVSc=HGVSc,
                            HGVSp=HGVSp,
                            allele_fraction=af,
                            mutation_type=consequence,
                            category=CATEGORY_SNV_INDEL,
                            gene=gene,
                        ).to_dict()
                    )

    return records


def find_drivers(args):
    hartwig = HartwigDecisionTree(args["gene_panel"], args["hotspot"])
    hartwig.build(args["inputvcf"], args["output"])


def dump_to_json(args):
    results = build_from_vcf(args["vcf"])
    cnvs = build_from_tsv(args["tsv_cnv"])
    results += cnvs
    json_object = json.dumps(results, indent=4)
    # Writing to sample.json
    with open(args["output"], "w") as outfile:
        outfile.write(json_object)


################################################
#   MAIN
################################################
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Tools for somatic variants annotation")
    subparsers = parser.add_subparsers(dest="command")

    driver_catalog = subparsers.add_parser("driverCatalogVCF", help="Create a catalog of drivers for SNVs and INDELs using the Hartwig Medical Foundation gene panel.")

    driver_catalog.add_argument("-i", "--inputvcf", help="somatic vcf containing normal and tumor genotypes. Order of the genotype columns: TUMOR, NORMAL.", required=True)
    driver_catalog.add_argument("-g", "--gene_panel", help="gene panel configuration", required=True)
    driver_catalog.add_argument(
        "-o",
        "--output",
        help="output VCF containing reported putative drivers",
        required=True,
    )

    driver_catalog.add_argument(
        "-s",
        "--hotspot",
        help="VCF containing hotspot mutations",
        required=True,
    )

    dump_json = subparsers.add_parser("dumpJSON", help="Dump somatic drivers into the JSON format.")

    dump_json.add_argument("-v", "--vcf", help="VCF file containing annotated putative drivers, must be annotated by the driverCatalogVCF command", required=True)
    dump_json.add_argument("-c", "--tsv_cnv", help="TSV file containing somatic CNVs. Should contain the following columns: chr, start, end, gene, category", required=True)
    dump_json.add_argument(
        "-o",
        "--output",
        help="output JSON containing reported putative drivers",
        required=True,
    )

    args = parser.parse_args()
    if args.command == "driverCatalogVCF":
        find_drivers(vars(args))

    if args.command == "dumpJSON":
        dump_to_json(vars(args))

import csv
import json
import argparse
from granite.lib.vcf_parser import Vcf
import pandas
import subprocess
import gzip


# VCF TAGS
#####################################################################################
CSQ_tag = "CSQ"
#####################################################################################

# VEP OUTPUT
# Definitions of the output fields available under "Default VEP output" section:
# https://www.ensembl.org/info/docs/tools/vep/vep_formats.html?redirect=no
#####################################################################################
FEATURE_field = "Feature"
FEATURE_TYPE_field = "Feature_type"
SYMBOL_field = "SYMBOL"
GENE_field = "Gene"
TRANSCRIPT = "Transcript"
CONSEQUENCE_field = "Consequence"
HGVSc_field = "HGVSc"
HGVSp_field = "HGVSp"
CANONICAL_field = "CANONICAL"
DRIVER_field = "DRIVER"
HOTSPOT_field = "HOTSPOT"
#####################################################################################

# OTHER CONSTANTS
#####################################################################################
CATEGORY_SNV_INDEL = "mutational"
CANONICAL_TRUE = "YES"
#####################################################################################

class Driver:
    """
    A class used to represent a driver mutation

    :param category : driver category
    :param gene: gene symbol

    """

    def __init__(self, gene, category=None):
        self.category = category
        self.gene = gene

    def to_dict(self):
        pass


class DriverSnvIndel(Driver):
    """
    A class used to represent a driver mutation for SNVs and INDELs

    :param vnt_obj : variant object
    :type vnt_obj: Variant

    :param transcript_consequence: transcript consequence
    :type transcript_consequence: str

    :param protein_mutation: protein mutation
    :type protein_mutation: str

    :param ens_gene: ensembl transcript gene ID
    :type ens_gene: str

    :param end_st: ensembl transcript protein ID
    :type end_st: str

    :param mutation_type: mutation type
    :type mutation_type: str

    :param hotspot: boolean value if it is a hotspot mutation
    :type hotspot: str

    """

    def __init__(
        self,
        vnt_obj,
        transcript_consequence=None,
        protein_mutation=None,
        ens_gene=None,
        transcript_id=None,
        mutation_type=None,
        hotspot=False,
        **kwargs,
    ):

        self.chrom = vnt_obj.CHROM
        self.pos = vnt_obj.POS
        self.ref = vnt_obj.REF
        self.alt = vnt_obj.ALT
        self.ens_gene = ens_gene
        self.transcript_id = transcript_id
        self.transcript_consequence = transcript_consequence
        self.protein_mutation = protein_mutation
        self.mutation_type = mutation_type
        self.hotspot = hotspot

        # tumor sample should be the first in the genotypes
        tumor_sample = vnt_obj.IDs_genotypes[0]
        self.allele_fraction = vnt_obj.get_genotype_value(tumor_sample, "AF")

        super().__init__(category="mutational", **kwargs)

    def to_dict(self):
        """

        A function used to serialize the object into a dictionary.
        It is used to dump the variant into a JSON format, some of the parameters are not included.
        All the None values are discarded.

        """
        fields = vars(self)
        return {
            k: v
            for (k, v) in fields.items()
            if k not in ["ens_gene", "transcript_id", "hotspot"] and v != None
        }

    def __str__(self):
        """

        String represation of the object, used for the DRIVER field in INFO.

        """
        return f"{self.gene}|{self.ens_gene}|{self.transcript_id}"


class DriverCnv(Driver):
    """
    A class used to represent a driver mutation for CNVs

    :param chr : chromosome
    :type chr: str

    :param start: start position of the segment
    :type start: str

    :param end: end position of the segment
    :type end: str

    """

    def __init__(self, chrom, start, end, **kwargs):
        self.chrom = chrom
        self.start = start
        self.end = end
        super().__init__(**kwargs)

    def to_dict(self):
        """

        A function used to serialize the object into a dictionary.
        It is used to dump the variant into the JSON format.

        """
        return vars(self)


class HartwigDecisionTree:
    """
    A class used to represent a Hartwig Decision Tree

    :param driver_panel: path to the driver panel used by HMF
    :type driver_panel: str

    :param hotspot_mutations: path to the VCF file containg hotspot mutations
    :type hotspot_mutations: str

    """

    def __init__(self, driver_panel, hotspot_mutations):

        # load the driver panel as a data frame
        gene_panel = pandas.read_csv(driver_panel, sep="\t")

        # create a list of genes for which missense variants should be reported
        self.report_missense = [
            gene[0]
            for gene in gene_panel.loc[
                gene_panel["reportMissense"] == True, ["gene"]
            ].values.tolist()
        ]

        # create a list of genes for which nonsense variants should be reported
        self.report_nonsense = [
            gene[0]
            for gene in gene_panel.loc[
                gene_panel["reportNonsense"] == True, ["gene"]
            ].values.tolist()
        ]

        # create a dictionary of genes for which additional reported transcripts should be reported
        # Keys: gene symbols
        # Values: lists of transcript IDs
        # Example:
        # {"CDKN2A" : ["ENST00000579755"] }
        self.non_canonical = (
            gene_panel.set_index("gene")["additionalReportedTranscripts"]
            .dropna()
            .to_dict()
        )

        # create a list of genes for which somatic hotspot mutations should be reported 
        self.hotspot = [
            gene[0]
            for gene in gene_panel.loc[
                gene_panel["reportSomaticHotspot"] == True, ["gene"]
            ].values.tolist()
        ]

        # create a dictionary of hotspot mutations
        # Keys: Variant key CHROM_POS_REF_ALT
        # Values: Variant object
        # Example:
        # {'chr16_2054395_CTC_AGG': granite.lib.vcf_parser.Vcf.Variant}
        self.hotspot_dict = {}
        hotspot_vcf = Vcf(hotspot_mutations)

        for vnt_obj in hotspot_vcf.parse_variants():
            self.hotspot_dict[self.__create_hotspot_key(vnt_obj)] = vnt_obj

    def __create_hotspot_key(self, vnt_obj):
        """
        Create an unique key for the variant: CHROM_POS_REF_ALT

        :param vnt_obj: variant
        :type vnt_obj: Variant

        :return: key for the variant
        :rtype: string
        """
        return f"{vnt_obj.CHROM}_{vnt_obj.POS}_{vnt_obj.REF}_{vnt_obj.ALT}"

    def build(self, input_vcf, save_vcf=None):
        """
        Find putative drivers mutation based on the HMF decision tree.
        Used for SNVs and INDELs.

        :param input_vcf: path to the input
        :type input_vcf: string

        :param save_vcf: path to the output VCF
        :type save_vcf: string

        :return: putative drivers
        :rtype: list
        """

        
        # Helper sub-function
        def check_non_canonical(gene, transcript):
            """
            Check if the given transcript should be reported for the given gene

            :param gene: gene symbol
            :type gene: str

            :param transcript: transcript ID
            :type transcript: str

            :return: True if the transcript should be reported, otherwise False
            :rtype: bool             
            """

            if gene in self.non_canonical.keys():
                if transcript in self.non_canonical[gene]:
                    return True
            return False

        # Helper sub-function   
        def check_hotspot(vnt_obj):
            """
            Check if the given variant is a hotspot mutation

            :param vnt_obj: variant to check
            :type vnt_obj: Variant

            :return: True if the variant is a hotspot mutation, otherwise False
            :rtype: bool
            """

            vnt_key = self.__create_hotspot_key(vnt_obj)
            if vnt_key in self.hotspot_dict.keys():
                return True
            else:
                return False

        # Helper sub-function
        def run_tree(vnt_obj):
            """
            Determine if a variant is a putative driver mutation

            :param vnt_obj: variant
            :type vnt_obj: Variant

            :return: DriverSnvIndel objects
            :rtype: list
            """

            # get a dictionary of transcripts from the CSQ tag
            transcripts_dict = vcf_obj.create_transcripts_dict(vnt_obj)

            # results list
            drivers = []

            #iterate over transcripts
            for transcript, values in transcripts_dict.items():

                ADD_TRANSCRIPT = False # flag to determine if the transcript should considered
                
                consequence = values[CONSEQUENCE_field] # Consequence value
                gene = values[SYMBOL_field] # SYMBOL value
                ens_gene = values[GENE_field] # Gene value
                HGVSc = values[HGVSc_field] # HGVSc value
                HGVSp = values[HGVSp_field] # HGVSp value
                transcript_id = values[FEATURE_field] # FEATURE value
                canonical = values[CANONICAL_field] # CANONICAL value
                
                
                # DECISION TREE LOGIC

                # consider all the canonical transcripts
                if canonical == CANONICAL_TRUE:
                    ADD_TRANSCRIPT = True
                else:
                    # check if non canonical transcript should be considered
                    ADD_TRANSCRIPT = check_non_canonical(gene, transcript)
                if (

                    # check if it is a missense variant and if the gene should be reported
                    (
                        ("missense_variant" in consequence)
                        and (gene in self.report_missense)
                    )
                    or (
                         # check if it is a nonsense variant and if the gene should be reported
                        (
                            "frameshift_variant" in consequence
                            or "stop_gained" in consequence
                        )
                        and (gene in self.report_nonsense)
                    )
                ) and ADD_TRANSCRIPT == True:
                    
                    # create a driver object and add to the results
                    drivers.append(
                        DriverSnvIndel(
                            vnt_obj=vnt_obj,
                            gene=gene,
                            ens_gene=ens_gene,
                            transcript_id=transcript_id,
                            transcript_consequence=HGVSc.split(":")[1],
                            protein_mutation=HGVSp.split(":")[1],
                        )
                    )

            # return list of putatvie drivers
            return drivers

        # Main
        vcf_obj = VcfVep(input_vcf)
        all_drivers = []

        # if return only driver objects
        if save_vcf == None:
            for vnt_obj in vcf_obj.parse_variants():
                all_drivers += run_tree(vnt_obj)

        # save to a VCF file
        else:

            with open(save_vcf, "w") as output:

                # adding header definitions

                # DRIVER definition
                vcf_obj.header.add_tag_definition(
                    f'##INFO=<ID={DRIVER_field},Number=.,Type=String,Description="Putative driver mutation calculated based on Hartwig\'s decision tree. Format: Gene|Ensembl Transcript ID|Protein change">',
                    tag_type="INFO",
                )
                
                # HOTSPOT definition
                vcf_obj.header.add_tag_definition(
                    f'##INFO=<ID={HOTSPOT_field},Number=.,Type=String,Description="Somatic hotspot location, contains a list of potentially affected genes">',
                    tag_type="INFO",
                )

                # writing the updated header to the file
                vcf_obj.write_header(output)

                # iterate over each variant in the file and find putative drivers
                for vnt_obj in vcf_obj.parse_variants():
                    # Find putative drivers
                    drivers = run_tree(vnt_obj)
                    all_drivers += drivers

                    # Add the DRIVER field for the variant only if there are drivers
                    if len(drivers) > 0:
                        DRIVER = ",".join([str(d) for d in drivers])
                        vnt_obj.add_tag_info(f"{DRIVER_field}={DRIVER}")

                    # Check if it is a hotspot mutation
                    hotspot = check_hotspot(vnt_obj)

                    if hotspot:

                        # list of gene symbols that should be reported as affected genes by the hotspot mutation
                        hotspot_genes = []

                        # get transcripts for the variant
                        transcripts_dict = vcf_obj.create_transcripts_dict(vnt_obj)

                        for _, values in transcripts_dict.items():

                            # check if a hotspot mutation should be reported for the given gene                      
                            if values[SYMBOL_field] in self.hotspot:

                                # add gene
                                hotspot_genes += [values[SYMBOL_field]]

                        #remove repetetive genes
                        hotspot_genes = set(hotspot_genes)

                        # Add the HOTSPOT field for the variant only if there are genes afftected by the hotspot mutation
                        if len(hotspot_genes) > 0:
                            vnt_obj.add_tag_info(
                                f"{HOTSPOT_field}={'|'.join(hotspot_genes)}"
                            )

                    # writing the varint to the VCF file
                    vcf_obj.write_variant(output, vnt_obj)

            # compression and indexing
            subprocess.run(["bgzip", save_vcf])
            subprocess.run(["tabix", save_vcf + ".gz"])

        return all_drivers


class VcfVep(Vcf):
    """
    An extension class of the Vcf class.
    This class is capable of creating a dictionary of transcripts for the given variant.
    """

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

    def create_transcripts_dict(self, vnt_obj):
        """
        Create a dictionary of transcripts for the variant.
        Each key is a transcript identifier, the value is a dictionary
        that holds all the information about the transcript.

        :param vnt_obj: variant
        :type vnt_obj: Variant

        :return: transcripts
        :rtype: dict
        """
        # get CSQ values
        fields = vnt_obj.get_tag_value(CSQ_tag).split(",")
        transcripts_dict = {}
        for field in fields:
            values = field.split("|")

            # merge CSQ fields with the values and convert it into a dictionary
            record = dict(zip(self.header.CSQ_fields, values))

            # add transcripts only to the final dictionary
            if record[FEATURE_TYPE_field] == TRANSCRIPT:
                transcript_id = record[FEATURE_field]
                transcripts_dict[transcript_id] = record
        return transcripts_dict

    class Header(Vcf.Header):
        def __init__(self, *args, **kwargs):
            super().__init__(*args, **kwargs)
            self.CSQ_fields = self.get_CSQ_fields()

        def get_CSQ_fields(self):
            """
            Get the CSQ fields from the header

            Example:
            For a VCF having the following definition of CSQ:
            ##INFO=<ID=CSQ,Number=.,Type=String,Description="Consequence annotations from Ensembl VEP. Format: Allele|Consequence|IMPACT|SYMBOL|Gene|Feature_type|Feature|BIOTYPE"

            This function will return: [Allele, Consequence, IMPACT, SYMBOL, Gene, Feature_type, Feature, BIOTYPE]

            :return: fields
            :rtype: list
            """
            for line in self.definitions.split("\n")[:-1]:
                # finding CSQ definition in the header
                if line.startswith("##INFO=<ID=CSQ,"):
                    try:
                        # getting format
                        format = line.split("Format:")[1]
                        # Cleaning format
                        format = (
                            format.replace("'", "").replace('"', "").replace(">", "")
                        )

                        # splitting fields
                        return format.split("|")
                    except Exception:
                        raise ValueError(
                            "Wrong format of the CSQ defitinion in the header."
                        )
            raise Exception("Missing CSQ definition in the header.")


def build_from_tsv(input_tsv):
    """
    Create a list of CNVs drivers from a TSV file.

    :param input_tsv: path to the input file
    :type input_tsv: str

    :return: drivers
    :rtype: list
    """
    # read a compressed TSV file 
    input_cnvs = csv.DictReader(gzip.open(input_tsv, mode="rt"), delimiter="\t")
    results = []
    # expected values in the TSV file
    fields = ["chrom", "start", "end", "gene", "category"]
    # iterate over each CNV
    for cnv in input_cnvs:
        
        # convert to a dictionary Example: { "chrom" : 1, "start": 100, "end ": 1000, "gene": TP53, "category": amplicitaion}
        record = dict(cnv)

        if sorted(record.keys()) != sorted(fields):
            raise Exception(
                f"Wrong fields in the CNV file expected {fields}, got {list(record.keys())}"
            )
        
        # cerate DriverCnv objectsd
        results += [DriverCnv(**record).to_dict()]
    return results


def build_from_vcf(input_vcf):
    """
    Create a list of SNVs and INDELs drivers from a VCF file.

    :param input_vcf: path to the input file
    :type input_vcf: str

    :return: drivers
    :rtype: list
    """
    records = []
    vcf_obj = VcfVep(input_vcf)

    # iterate over variants
    for vnt_obj in vcf_obj.parse_variants():

        # get a dictionary of transcripts from the CSQ tag
        transcripts_dict = vcf_obj.create_transcripts_dict(vnt_obj)

        # list of genes affected by hotspot mutations
        hotspot_genes = []
        # try to get HOTSPOT from INFO
        try:
            hotspot_genes = vnt_obj.get_tag_value(HOTSPOT_field).split("|")

            # iterate over transcripts, find canonical one, and collect all the necessary information
            for _, values in transcripts_dict.items():

                # check if for the given gene hotspot mutation should be reported and check if it is a canonical transcript
                if (
                    values[SYMBOL_field] in hotspot_genes
                    and values[CANONICAL_field] == CANONICAL_TRUE
                ):
                    # create a driver object, add to the final list
                    records.append(
                        DriverSnvIndel(
                            vnt_obj=vnt_obj,
                            ens_gene=values[GENE_field],
                            transcript_id=values[FEATURE_field],
                            transcript_consequence=values[HGVSc_field].split(":")[1],
                            protein_mutation=values[HGVSp_field].split(":")[1],
                            mutation_type=values[CONSEQUENCE_field],
                            gene=values[SYMBOL_field],
                        ).to_dict()
                    )
        # TODO: we should change it into custom errors
        
        except Exception:
            # pass if it is not a hotspot mutation
            pass
        
        # try to get DRIVER from INFO
        try:
            drivers = vnt_obj.get_tag_value(DRIVER_field).split(",")

        # TODO: we should change it into custom error
        except Exception:
            # continue the loop, go to the next variant if DRIVER is missing
            continue

        for driver in drivers:
            gene_d, _, transcript_id_d = driver.split("|")
            # if its in hotspot_genes, it's been already added
            if gene_d not in hotspot_genes:
                # get values for the transcript ID
                trans = transcripts_dict[transcript_id_d]

                # create a driver object, add to the final list
                records.append(
                    DriverSnvIndel(
                        vnt_obj=vnt_obj,
                        ens_gene=trans[GENE_field],
                        transcript_id=trans[FEATURE_field],
                        transcript_consequence=trans[HGVSc_field].split(":")[1],
                        protein_mutation=trans[HGVSp_field].split(":")[1],
                        mutation_type=trans[CONSEQUENCE_field],
                        gene=trans[SYMBOL_field],
                    ).to_dict()
                )

    # return all the driver objects
    return records


def find_drivers(args):
    """
    Run Hartwig Decision Tree to find putative drivers
    and save to a new VCF file

    :param args: arguments from the CLI
    :type args: dict
    """
    hartwig = HartwigDecisionTree(args["gene_panel"], args["hotspot"])
    hartwig.build(args["inputvcf"], args["output"])


def dump_to_json(args):
    """
    Dump driver mutations (SNVs, INDELs, CNVs) into a JSON file.

    :param args: arguments from the CLI
    :type args: dict
    """
    results = build_from_vcf(args["vcf"])
    cnvs = build_from_tsv(args["tsv_cnv"])
    results += cnvs
    json_object = json.dumps(results, indent=4)
    with open(args["output"], "w") as outfile:
        outfile.write(json_object)


################################################
#   MAIN
################################################
if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Tools for somatic variants annotation"
    )
    subparsers = parser.add_subparsers(dest="command")

    driver_catalog = subparsers.add_parser(
        "driverCatalogVCF",
        help="Create a catalog of drivers for SNVs and INDELs using the gene panel used by Hartwig Medical Foundation.",
    )

    driver_catalog.add_argument(
        "-i",
        "--inputvcf",
        help="somatic vcf containing normal and tumor genotypes. Order of the genotype columns: TUMOR, NORMAL.",
        required=True,
    )
    driver_catalog.add_argument(
        "-g", "--gene_panel", help="gene panel configuration", required=True
    )
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

    dump_json = subparsers.add_parser(
        "dumpJSON", help="Dump somatic drivers into the JSON format."
    )

    dump_json.add_argument(
        "-v",
        "--vcf",
        help="VCF file containing annotated putative drivers, must be annotated by the driverCatalogVCF command",
        required=True,
    )
    dump_json.add_argument(
        "-c",
        "--tsv_cnv",
        help="TSV file containing somatic CNVs. Should contain the following columns: chr, start, end, gene, category",
        required=True,
    )
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

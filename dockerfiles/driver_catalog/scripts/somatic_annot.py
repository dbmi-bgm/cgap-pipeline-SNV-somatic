import csv
import json
import argparse
from granite.lib.vcf_parser import Vcf
import pandas
import subprocess
import gzip

CSQ_tag = "CSQ"
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
CATEGORY_SNV_INDEL = "mutational"
CANONICAL_TRUE = "YES"


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

        tumor_sample = vnt_obj.IDs_genotypes[0]
        self.allele_fraction = vnt_obj.get_genotype_value(tumor_sample, "AF")

        super().__init__(category="mutational", **kwargs)

    def to_dict(self):
        """

        A function used to serialize the object into a dictionary.
        It is used to dump the varian into a JSON format, some of the parameters are not included.
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
        It is used to dump the varian into a JSON format.

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

        self.non_canonical = (
            gene_panel.set_index("gene")["additionalReportedTranscripts"]
            .dropna()
            .to_dict()
        )
        self.hotspot = [
            gene[0]
            for gene in gene_panel.loc[
                gene_panel["reportSomaticHotspot"] == True, ["gene"]
            ].values.tolist()
        ]
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

        # Helper functions

        def check_non_canonical(gene, transcript):

            if gene in self.non_canonical.keys():
                if transcript in self.non_canonical[gene]:
                    return True
            return False

        def query_hotspot(vnt_obj):
            vnt_key = self.__create_hotspot_key(vnt_obj)
            if vnt_key in self.hotspot_dict.keys():
                return True
            else:
                return False

        def query(vnt_obj):
            transcripts = vnt_obj.get_tag_value(CSQ_tag).split(",")
            transcripts_dict = vcf_obj.create_transcripts_dict(vnt_obj)
            drivers = []

            for transcript, values in transcripts_dict.items():
                ADD_TRANSCRIPT = False
                consequence = values[CONSEQUENCE_field]
                gene = values[SYMBOL_field]
                ens_gene = values[GENE_field]
                HGVSc = values[HGVSc_field]
                HGVSp = values[HGVSp_field]
                transcript_id = values[FEATURE_field]
                canonical = values[CANONICAL_field]

                if canonical == CANONICAL_TRUE:
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
            return drivers

        # Main
        vcf_obj = VcfVep(input_vcf)
        all_drivers = []
        if save_vcf == None:
            for vnt_obj in vcf_obj.parse_variants():
                all_drivers += query(vnt_obj)

        else:
            with open(save_vcf, "w") as output:
                vcf_obj.header.add_tag_definition(
                    f'##INFO=<ID={DRIVER_field},Number=1,Type=String,Description="Putative driver mutation calculated based on Hartwig\'s decision tree. Format: Gene|Ensemble Transcript ID|Protein change">',
                    tag_type="INFO",
                )
                vcf_obj.header.add_tag_definition(
                    f'##INFO=<ID={HOTSPOT_field},Number=1,Type=String,Description="Somatic hotspot location. Format: Gene|Ensemble Transcript ID|Protein change ">',
                    tag_type="INFO",
                )
                vcf_obj.write_header(output)
                for vnt_obj in vcf_obj.parse_variants():
                    # Find putative drivers
                    drivers = query(vnt_obj)
                    all_drivers += drivers

                    # Add the DRIVER field for the variant only if there are drivers
                    if len(drivers) > 0:
                        DRIVER = ",".join([str(d) for d in drivers])
                        vnt_obj.add_tag_info(f"{DRIVER_field}={DRIVER}")

                    # Check if it is a hotspot mutation
                    hotspot = query_hotspot(vnt_obj)

                    if hotspot:
                        hotspot_genes = []
                        transcripts_dict = vcf_obj.create_transcripts_dict(vnt_obj)

                        for _, v in transcripts_dict.items():
                            if (
                                v[SYMBOL_field] in self.hotspot
                                and v[SYMBOL_field] not in hotspot_genes
                            ):
                                hotspot_genes += [v[SYMBOL_field]]

                        if len(hotspot_genes) > 0:
                            vnt_obj.add_tag_info(
                                f"{HOTSPOT_field}={'|'.join(hotspot_genes)}"
                            )

                    vcf_obj.write_variant(output, vnt_obj)

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
        transcripts = vnt_obj.get_tag_value(CSQ_tag).split(",")
        transcripts_dict = {}
        for transcript in transcripts:
            transcript_fields = transcript.split("|")

            record = dict(zip(self.header.CSQ_fields, transcript_fields))
            if record[FEATURE_TYPE_field] == TRANSCRIPT:
                transcript_id = record[FEATURE_field]
                transcripts_dict[transcript_id] = record
        return transcripts_dict

    class Header(Vcf.Header):
        def __init__(self, *args, **kwargs):
            super().__init__(*args, **kwargs)
            self.CSQ_fields = self.get_CSQ_fields()

        def get_CSQ_fields(self):
            for line in self.definitions.split("\n")[:-1]:
                if line.startswith("##INFO=<ID=CSQ,"):
                    try:
                        format = line.split("Format:")[1]
                        # Cleaning format
                        format = (
                            format.replace("'", "").replace('"', "").replace(">", "")
                        )

                        return format.split("|")
                    except Exception:
                        raise ValueError(
                            "Wrong format of the CSQ defitinion in the header."
                        )


def build_from_tsv(input_tsv):
    """
    Create a list of CNVs drivers from a TSV file.

    :param input_tsv: path to the input file
    :type input_tsv: str

    :return: drivers
    :rtype: list
    """
    input_file = csv.DictReader(gzip.open(input_tsv, mode="rt"), delimiter="\t")
    results = []
    fields = ["chrom", "start", "end", "gene", "category"]
    for inn in input_file:
        record = dict(inn)
        if sorted(record.keys()) != sorted(fields):
            raise Exception(
                f"Wrong fields in the CNV file expected {fields}, got {list(record.keys())}"
            )
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
    print(vcf_obj.header.check_tag_definition("CSQ"))
    for vnt_obj in vcf_obj.parse_variants():
        transcripts_dict = vcf_obj.create_transcripts_dict(vnt_obj)

        hotspot_genes = []
        try:
            hotspot_genes = vnt_obj.get_tag_value(HOTSPOT_field).split("|")
            for k, v in transcripts_dict.items():
                if v[SYMBOL_field] in hotspot_genes and v[CANONICAL_field] == CANONICAL_TRUE:
                    records.append(
                        DriverSnvIndel(
                            vnt_obj=vnt_obj,
                            ens_gene=ens_gene,
                            transcript_id=v[FEATURE_field],
                            transcript_consequence=v[HGVSc_field].split(":")[1],
                            protein_mutation=v[HGVSp_field].split(":")[1],
                            mutation_type=v[CONSEQUENCE_field],
                            gene=v[SYMBOL_field],
                        ).to_dict()
                    )
        except Exception:
            pass

        try:
            drivers = vnt_obj.get_tag_value(DRIVER_field).split(",")
        except Exception:
            continue

        for driver in drivers:
            gene_d, transcript_gene_d, transcript_id_d = driver.split("|")
            # if its in hotspot_genes, it's been already added
            if gene_d not in hotspot_genes:
                trans = transcripts_dict[transcript_id_d]
                records.append(
                    DriverSnvIndel(
                        vnt_obj=vnt_obj,
                        ens_gene=ens_gene,
                        transcript_id=trans[FEATURE_field],
                        transcript_consequence=trans[HGVSc_field].split(":")[1],
                        protein_mutation=trans[HGVSp_field].split(":")[1],
                        mutation_type=trans[CONSEQUENCE_field],
                        gene=trans[SYMBOL_field],
                    ).to_dict()
                )

    return records


def find_drivers(args):
    """
    Run Hartwig Decition Tree to find putative drivers
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

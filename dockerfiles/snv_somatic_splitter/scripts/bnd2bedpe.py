import argparse
import csv
from granite.lib import vcf_parser

SVTYPE = "SVTYPE"
MATEID = "MATEID"
CIPOS = "CIPOS"
BND = "BND"


INFO_fields = ["IMPRECISE"]


class BedpeRecord:

    def __init__(self, chrom1="", start1="", end1="", chrom2="", start2="", sv_id="", pe_support="",
                 strand1="", strand2="", svclass=""):
        """
        BedpeRecord represents a single bedpe record
        @param chrom1: chrom1 field in the BEDPE format
        @param start1: start1 field in the BEDPE format
        @param end1:  end1 field in the BEDPE format
        @param chrom2: chrom2 field in the BEDPE format
        @param start2: start2 field in the BEDPE format
        @param sv_id: id of structural variant
        @param pe_support: phred-scaled quality score
        @param strand1: strand for the first end of the feature
        @param strand2: strand for the second end of the feature
        @param svclass: structural variant type
        """
        self.chrom1 = chrom1
        self.start1 = start1
        self.end1 = end1
        self.chrom2 = chrom2
        self.start2 = start2
        self.sv_id = sv_id
        self.pe_support = pe_support
        self.strand1 = strand1
        self.strand2 = strand2
        self.svclass = svclass


class Bnd:
    def __init__(self, vnt_obj):
        """
        Bnd represents a single breakend
        @rtype: object
        @param vnt_obj: Variant class object
        """
        self.CHROM = vnt_obj.CHROM
        self.POS = vnt_obj.POS
        self.ID = vnt_obj.ID
        self.REF = vnt_obj.REF
        self.ALT = vnt_obj.ALT.split(",")
        self.QUAL = vnt_obj.QUAL
        self.INFO = vnt_obj.INFO
        self.MATES = self.get_mates(vnt_obj)
        self.SGL = False  # single breakend

        if len(self.MATES) > 0:
            strands = [
                "-" if any(alt.startswith(x) for x in ["[", "]"]) else "+"
                for alt in vnt_obj.ALT.split(",")
            ]
            self.STRANDS = dict(zip(self.MATES, strands))
        else:
            self.SGL = True
            if self.ALT[0].startswith("."):  # TODO: check how many alternates
                self.STRANDS = {None: "-"}

            elif self.ALT[0].endswith("."):
                self.STRANDS = {None: "+"}

        self.CIPOS = self._get_cipos(vnt_obj)

    def calculate_bedpe_end(self):
        """
        Calculate the one-based ending position of the end of the feature
        @return: Position
        """
        if self.CIPOS:
            return self.POS + abs(self.CIPOS[1])
        else:
            return self.POS

    def calculate_bedpe_start(self):
        """
        Calculate the zero-based starting position of the end of the feature
        @return: Position
        """
        if self.CIPOS:
            return self.POS - self.CIPOS[0] - 1
        else:
            return self.POS - 1

    def get_mates(self, vnt_obj):
        """
        List all mates of the variant
        @param vnt_obj: Variant object
        @return: List of the mates. If there are no mates, returns an empty list.
        """
        try:
            return vnt_obj.get_tag_value(MATEID).split(",")
        except:
            return []

    def _get_cipos(self, vnt_obj):
        """
        Get values of the CIPOS field
        @param vnt_obj: Variant objectt
        @return: List of the CIPOS values
        """
        try:
            return list(map(int, vnt_obj.get_tag_value(CIPOS).split(",")))
        except:
            return None


class Converter:
    def __init__(self, bnds):
        """
        Converter class is responsible for the conversion of SVs from VCF to BEDPE
        @param bnds: a list of variant ID pairs that are mates
        """
        self.bnd_ids = bnds.keys()
        self.bnds = bnds
        self.pairs = ()

    def pair_bnds(self):
        """
        Create a list of variant ID pairs that are mates
        @rtype: object
        """
        pairs = []
        for _, bnd_obj in self.bnds.items():
            pairs += [
                [bnd_obj.ID, mate_id]
                for mate_id in bnd_obj.MATES
                if ([mate_id, bnd_obj.ID] not in pairs)
            ]
            if len(bnd_obj.MATES) == 0:
                pairs += [[bnd_obj.ID, None]]


        self.pairs = pairs
        return self.pairs

    def create_bedpe_records(self):
        """
        Convert variants into the BEDPE format
        @return:
        """
        records = []
        chromosomes = [str(chrom) for chrom in list(range(1, 23))] + ["X", "Y"]
        chromosomes += ["chr" + chrom for chrom in chromosomes]
        for pair in self.pairs:
            if not self.bnds[pair[0]].SGL:
                record = self.create_paired_record(pair)

            else:
                record = self.create_single_breakend_record(pair)
            records.append(record)

        return records

    def create_single_breakend_record(self, pair):
        """
        Convert a single breakend record into a dictionary with keys that correspond to the BEDPE format
        @param pair:
        @return: BEDPE record as dictionary
        """
        breakend = self.bnds[pair[0]]
        record = BedpeRecord()
        record.chrom1 = breakend.CHROM
        record.start1 = int(breakend.POS) - 1
        record.end1 = breakend.calculate_bedpe_end()
        record.chrom2 = "."
        record.start2 = "-1"
        record.end2 = "-1"
        record.sv_id = breakend.ID
        record.pe_support = str((float(breakend.QUAL)))

        record.strand1 = breakend.STRANDS[None]
        record.strand2 = "."
        record.svclass = "."

        return record

    def create_paired_record(self, pair):
        """
        Convert mated SVs into a dictionary with keys that correspond to the BEDPE format

        @param pair: Pair of variants IDs that are mates
        @return:
        """
        record = {}
        pair1 = self.bnds[pair[0]]
        pair2 = self.bnds[pair[1]]
        chrom1 = pair1.CHROM
        chrom2 = pair2.CHROM
        strand1 = pair1.STRANDS[pair[1]]
        strand2 = pair2.STRANDS[pair[0]]
        record.chrom1 = chrom1
        record.start1 = pair1.calculate_bedpe_start()
        record.end1 = pair1.calculate_bedpe_end()
        record.chrom2 = chrom2
        record.start2 = pair2.calculate_bedpe_start()
        record.end2 = pair2.calculate_bedpe_end()
        record.sv_id = pair1.ID
        record.pe_support = str((float(pair1.QUAL)))
        record.strand1 = strand1
        record.strand2 = strand2
        if chrom1 == chrom2:
            if strand1 == "+":
                if strand2 == "-":
                    svtype = "DEL"  # DEL = +-
                else:
                    svtype = "h2hINV"  # h2hINV = ++
            else:
                if strand2 == "+":
                    svtype = "DUP"  # DUP = -+
                else:
                    svtype = "t2tINV"  # t2tINV = --
        else:
            svtype = "TRA"

        record.svclass = svtype
        return record


def vcf_scanner(vcf_input):
    """
    Scans a VCF file and converts BND SVs into Bnd objects
    @param vcf_input: Path to the input VCF file
    @return: Dictionary of Bnd objects. Keys are IDs of the variants, values are Bnd objects.
    """
    vcf = vcf_parser.Vcf(vcf_input)
    bnd_dict = {}
    for vnt_obj in vcf.parse_variants():
        if vnt_obj.get_tag_value(SVTYPE) == BND:
            bnd_obj = Bnd(vnt_obj)
            bnd_dict[bnd_obj.ID] = bnd_obj
    return bnd_dict


def bnd2bedpe(args):
    bedpe_columns = [
        "chrom1",
        "start1",
        "end1",
        "chrom2",
        "start2",
        "end2",
        "sv_id",
        "pe_support",
        "strand1",
        "strand2",
        "svclass"
    ]
    bnd_dict = vcf_scanner(args["input"])
    converter = Converter(bnd_dict)
    converter.pair_bnds()
    records = converter.create_bedpe_records()
    csv_file = args["output"]
    with open(csv_file, "w") as csvfile:
        writer = csv.DictWriter(csvfile, fieldnames=bedpe_columns, delimiter="\t")
        writer.writeheader()
        for data in records:
            writer.writerow(vars(data))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="")
    parser.add_argument("-i", "--input", help="input SV VCF", required=True)
    parser.add_argument("-o", "--output", help="output SV BEDPE", required=True)
    args = vars(parser.parse_args())
    bnd2bedpe(args)

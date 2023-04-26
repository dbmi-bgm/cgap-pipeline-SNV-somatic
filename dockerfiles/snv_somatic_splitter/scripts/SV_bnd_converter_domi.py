#!/usr/bin/env python3

################################################
#
#  Script to convert Sentieon SVs to BEDPE
#
################################################


################################################
#   Libraries
################################################
import argparse
from granite.lib import vcf_parser
import csv


################################################
#   Classes
################################################
class Bnd:
    def __init__(self, vnt_obj, index):
        self.CHROM = vnt_obj.CHROM
        self.POS = vnt_obj.POS
        self.ID = vnt_obj.ID
        self.REF = vnt_obj.REF
        self.ALT = vnt_obj.ALT.split(",")
        self.QUAL = vnt_obj.QUAL
        self.INFO = vnt_obj.INFO
        self.MATES = self._mate_lister(vnt_obj)
        self.INDEX = index
        self.sgl = False  # single breakend
        if len(self.MATES) > 0:
            strands = [
                "-" if any(alt.startswith(x) for x in ["[", "]"]) else "+"
                for alt in vnt_obj.ALT.split(",")
            ]
            self.STRANDS = dict(zip(self.MATES, strands))
        else:
            self.sgl = True
            if self.ALT[0].startswith("."):  # TODO: check how many alternates
                self.STRANDS = {None: "-"}
            elif self.ALT[0].endswith("."):
                self.STRANDS = {None: "+"}

        """if len(self.MATES) == 0:
            print(self.REF)
            print(self.STRANDS)
            print(self.ALT)
            print("*************************")"""

    def _mate_lister(self, vnt_obj):
        try:
            return vnt_obj.get_tag_value("MATEID").split(",")
        except:
            return []


class Converter:
    def __init__(self, bnds):
        self.bnd_ids = bnds.keys()
        self.bnds = bnds
        self.pairs = ()

    def pair_bnds(self):
        pairs = []

        for _, bnd_obj in self.bnds.items():
            print(bnd_obj.CHROM)
            pairs += [
                [bnd_obj.ID, mate_id]
                for mate_id in bnd_obj.MATES
                if ([mate_id, bnd_obj.ID] not in pairs)
            ]

        self.pairs = pairs
        return self.pairs

    def create_bedpe_records(self):
        records = []
        chromosomes = [str(chrom) for chrom in list(range(1, 23))] + ["X", "Y"]
        chromosomes += ["chr" + chrom for chrom in chromosomes]

        for pair in self.pairs:
            pair1 = self.bnds[pair[0]]
            pair2 = self.bnds[pair[1]]

            # only want the main chromosomes
            if pair1.CHROM not in chromosomes or pair2.CHROM not in chromosomes:
                pass
            else:
                record = {}
                CHROM1 = pair1.CHROM
                CHROM2 = pair2.CHROM
                STRAND1 = pair1.STRANDS[pair[1]]
                STRAND2 = pair2.STRANDS[pair[0]]
                record["chrom1"] = CHROM1
                record["start1"] = str(int(pair1.POS) - 1)
                record["end1"] = str(pair1.POS)
                record["chrom2"] = CHROM2
                record["start2"] = str(int(pair2.POS) - 1)
                record["end2"] = str(pair2.POS)
                record["sv_id"] = pair1.ID
                record["pe_support"] = str(pair1.QUAL)
                record["strand1"] = STRAND1
                record["strand2"] = STRAND2

                if CHROM1 == CHROM2:
                    if STRAND1 == "+":
                        if STRAND2 == "-":
                            SVTYPE = "DEL"  # DEL = +-
                        else:
                            SVTYPE = "h2hINV"  # h2hINV = ++

                    else:
                        if STRAND2 == "+":
                            SVTYPE = "DUP"  # DUP = -+
                        else:
                            SVTYPE = "t2tINV"  # t2tINV = --
                else:
                    SVTYPE = "TRA"

                record["svclass"] = SVTYPE
                record["svmethod"] = "Sentieon_TNscope"

                records.append(record)

        return records


def vcf_scanner(vcf_input):
    vcf = vcf_parser.Vcf(vcf_input)
    bnd_dict = {}

    for index, vnt_obj in enumerate(vcf.parse_variants()):
        if vnt_obj.get_tag_value("SVTYPE") == "BND":
            bnd_obj = Bnd(vnt_obj, index)
            bnd_dict[bnd_obj.ID] = bnd_obj
        else:
            # so far, I have only seen SVTYPE=INS in tumor only analyses
            # we will only support paired breakpoints at this point with BEDPE
            if vnt_obj.get_tag_value("SVTYPE") == "INS":
                pass

    return bnd_dict


def main(args):
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
        "svclass",
        "svmethod",
    ]
    bnd_dict = vcf_scanner(args["inputvcf"])
    converter = Converter(bnd_dict)
    converter.pair_bnds()
    records = converter.create_bedpe_records()

    csv_file = args["outputBEDPE"]
    with open(csv_file, "w") as csvfile:
        writer = csv.DictWriter(csvfile, fieldnames=bedpe_columns, delimiter="\t")
        writer.writeheader()
        for data in records:
            writer.writerow(data)


################################################
#   MAIN
################################################
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="")

    parser.add_argument("-i", "--inputvcf", help="input SV vcf", required=True)
    parser.add_argument("-o", "--outputBEDPE", help="output bedPE", required=True)

    args = vars(parser.parse_args())

    main(args)

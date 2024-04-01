import pyranges
import argparse
import subprocess
from collections import defaultdict
import os

BEDTOOLS = "bedtools"

# 将bed文件拆分成仅包含外显子区域

# 只测试过人的数据!!

class ExonIntron:
    """
    """
    def __init__(self, gtf_file, out_dir, bed_file, keep_chr = None):
        if not os.path.exists(f"{out_dir}"):
            cmd = f'mkdir -p {out_dir}'
            subprocess.check_call(cmd, shell=True)
        self.keep_chr = keep_chr
        if keep_chr is None:
            self.keep_chr = "all" #[str(chrs + 1) for chrs in range(22)] + ["X", "Y"]

        self.gtf_file = gtf_file
        self.out_dir = out_dir
        self.bed_file = bed_file
        self.bed_out_name = self.bed_file.split('/')[-1].split(".")[0]
        self.exon_out = f"{self.out_dir}/exon_out.bed"
        self.intron_out = f"{self.out_dir}/{self.bed_out_name}_intron_out.bed"
        self.bed_out = f"{self.out_dir}/{self.bed_out_name}_without_intron.bed"


    def read_gtf(self):
        self.gtf_df = pyranges.read_gtf(self.gtf_file).as_df()
        # start + 1 后续bedtools无需-1
        self.gtf_df["Start"] = self.gtf_df["Start"] + 1
    
    def set_keep_chrs(self):
        if self.keep_chr == "all":
            self.keep_chr = list(set(self.gtf_df["Chromosome"].to_list()))
            self.keep_chr = list(map(lambda x: str(x), self.keep_chr))

        else:
            self.keep_chr = self.keep_chr.split(":")
            self.gtf_df = self.gtf_df[self.gtf_df["Chromosome"].isin(self.keep_chr)]

    def extract_exon(self):
        transcript_exon_dict = defaultdict(list)
        for _, row in self.gtf_df.iterrows():
            chrs = row["Chromosome"]
            entry_type = row["Feature"]
            transcript_id = row["transcript_id"]

            items = [
                str(row[attr]) for attr in [
                    "Strand", "Feature", "gene_id", "gene_name",
                    "gene_biotype", "transcript_id"
                ]
            ]

            record_k = ":".join(items)

            if entry_type == "exon" and str(chrs) in self.keep_chr:
                exon = (str(chrs) ,row["Start"], row["End"])
                if transcript_id not in transcript_exon_dict.keys():
                    transcript_exon_dict[record_k].append(exon)
            else:
                continue
        # sort
        for key, _ in transcript_exon_dict.items():
            transcript_exon_dict[key] = sorted(transcript_exon_dict[key], key = lambda x: x[1])

        self.transcript_exon_dict = transcript_exon_dict

    def extract_intron(self):
        region_intron_dict = defaultdict(list)

        for k, region_all in self.transcript_exon_dict.items():
            region_intron = []
            for n in range(len(region_all) - 1):
                intron = (region_all[n][0] ,region_all[n][2] + 1, region_all[n + 1][1])
                region_intron.append(intron)
            region_intron_dict[k] = region_intron

        self.region_intron_dict = region_intron_dict

    def write_to_bed(self, feature):
        if feature == "exon":
            fd = open(f"{self.exon_out}", "w")
            target_dict = self.transcript_exon_dict
        elif feature == "intron":
            fd = open(f"{self.intron_out}", "w")
            target_dict = self.region_intron_dict

        for info, region_all in target_dict.items():
            gene_info_lst = info.split(':')
            if "-" in gene_info_lst:
                region_all.sort(key = lambda x: x[1], reverse = True)

            for region in region_all:
                region = list(map(lambda x:str(x), region))
                record = "\t".join(region) + "\t" + "\t".join(gene_info_lst)
                fd.write(record)
                fd.write("\n")

    def remove_intron_from_bed(self):
        cmd = f"{BEDTOOLS} subtract -a {self.bed_file} -b {self.intron_out} > {self.bed_out}"
        subprocess.check_call(cmd, shell = True)


    def run(self):
        self.read_gtf()
        self.set_keep_chrs()
        self.extract_exon()
        self.extract_intron()
        self.write_to_bed("exon")
        self.write_to_bed("intron")
        self.remove_intron_from_bed()

def main():
    parser = argparse.ArgumentParser(description='')
    parser.add_argument(f'--gtf',help='', required=True)
    parser.add_argument(f'--bed',help='', required=True)
    parser.add_argument(f'--out_dir',help='', required=True)
    parser.add_argument(f'--keep_chrs',help='chrs1:chrs2:chrs3', required=False)
    args = parser.parse_args()
    target = ExonIntron(
        gtf_file = args.gtf,
        out_dir = args.out_dir,
        bed_file = args.bed,
        keep_chr = args.keep_chrs
        )
    target.run()

if __name__ == "__main__":
    main()
from collections import defaultdict
import numpy as np
import pandas as pd
import glob

def fa_iter(in_fa):
    with open(in_fa) as fd:
        while True:
            line = fd.readline()
            line = line.strip("\n")
            if not line:
                break
            yield line

def format(in_fa, out_fa):
    cache = 500000
    cache_dict = defaultdict(lambda:defaultdict(str))
    with open(out_fa, "w") as fd:
        for line in fa_iter(in_fa):
            if line.startswith(">"):
                seq_name = line.strip("\n")
                cache_dict.update({seq_name:[]})
            else:
                cache_dict[seq_name].append(line.strip("\n"))
            if len(cache_dict) == cache:
                for s_name, s in cache_dict.items():
                    seq = "".join(s)
                    fd.write(f"{s_name}\n{seq}\n")
        
        if len(cache_dict) != 0:
            for s_name, s in cache_dict.items():
                seq = "".join(s)
                fd.write(f"{s_name}\n{seq}\n")

def mk_region(fa):
    names_lst = []
    lens = []
    for line in fa_iter(fa):
        if line.startswith(">"):
            names_lst.append(line)
        else:
            lens.append(len(line))
    
    ends = np.cumsum(lens)
    inx = 1
    new_com = []
    for i, _ in enumerate(ends):
        if inx == 1:
            new_com.append((1, ends[i]))
            inx+=1
        else:
            new_com.append((ends[i-1]+1, ends[i]))
    
    # make region dict
    region_dict = {}
    ranges = 25
    for k, v in zip(names_lst, new_com):
        region_dict.update({k.replace('>', ''):v})
        for sufx in range(ranges):
            region_dict.update({f"{k.replace('>', '')}_{sufx}":v})
    return region_dict

def fix_loc(df, region_dict):
    records = [] # name start end
    for _, row in df.iterrows():
        new_start = row['start'] + region_dict[ row['GeneID']][0] 
        new_end = row['end'] + region_dict[ row['GeneID']][0]
        records.append((row['GeneID'], new_start, new_end))
    
    df = pd.DataFrame(records)
    df.columns = ["gene", "start", "end"]
    df = df.sort_values(by="start")
    return df

def main():
    fa_in = 'MAGs.S9.Bin76.fa' # 输入的fasta
    fa_out = "format.fa"       # 删除不必要换行符的fasta
    data_dir = "./test/"       # 注释文件目录

    files = glob.glob(f"{data_dir}/*.xls") # 默认xls结尾, 可以根据实际自己改
    format(fa_in, fa_out)
    region_dict = mk_region(fa_out)
    df_in = pd.concat([pd.read_table(f) for f in files ])
    df_out = fix_loc(df_in, region_dict)
    df_out.to_csv("./out_table.tsv", sep="\t", index=None)

if __name__ == '__main__':
    main()
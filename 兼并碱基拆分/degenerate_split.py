from collections import Counter
from tk.transfer import fasta_to_dict
import argparse


# 快速还原兼并碱基

# --fa 含有兼并碱基的fasta文件

DICT = {
    "H": ["A", "C", "T"],
    "B": ["G", "T", "C"],
    "R": ["A", "G"],
    "Y": ["C", "T"],
    "M": ["A", "C"],
    "K": ["T", "G"],
    "S": ["G", "C"],
    "W": ["A", "T"],
    "B": ["G", "T", "C"],
    "V": ["G", "A", "C"],
    "D": ["G", "A", "T"]
}

BASE = list(DICT.keys())


def get_dbase_list(seq):
    # 获取某条序列中的兼并碱基信息
    # 例如[H,H,B]
    seq_dict = Counter(seq)
    out = []
    for key, value in seq_dict.items():
        if key not in ["A", "T", "G", "C", "N", '-']:
            out += list(key * value)
    return out


def check_is_have_target(seqs: list, target):
    # 查看序列中是否还有兼并碱基
    flags = []
    for seq in seqs:
        for t in target:
            flags.append(seq.find(t))
    flag = list(set(flags)) != [-1]
    return flag

def _replace_seq(seqs, base):
    # 替换兼并碱基
    # 输入:
    # seqs: 序列列表
    # base: get_dbase_list的输出
    # 输出: 替换后的list
    # 每次替换一个兼并碱基
    r_bases = DICT[base]
    outs = []
    for seq in seqs:
        for r_base in r_bases:
            outs.append(seq.replace(base, r_base, 1))
    return outs


def seq_explosion(seqs: list):
    # 递归替换所有的兼并碱基

    if not check_is_have_target(seqs, BASE):
        return seqs
    else:
        new_seqs = _replace_seq(seqs, dbases.pop())
        print(new_seqs)
        return seq_explosion(new_seqs)

def main():
    parser = argparse.ArgumentParser(description='')
    parser.add_argument('--fa',"-f", help='', required=True)
    args = parser.parse_args()
    outs = []
    fa_dict = fasta_to_dict(args.fa)
    for seq, name in fa_dict.items():
        
        #  必须生成全局变量dbases
        global dbases
        
        dbases = get_dbase_list(seq)
        
        if len(dbases) == 0:
            outs.append(f">{name}\n{seq}\n")
        else:
            dseqs = seq_explosion([seq,])
            for inx, dseq in enumerate(dseqs):
                new_name = f"{name}_{inx+1}"
                outs.append(f">{new_name}\n{dseq}\n")
    # write
    with open("./output2.fa", "w") as fd:
        for r in outs:
            fd.write(r)

if __name__ == '__main__':
    main()
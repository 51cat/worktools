import glob
from scipy.io import mmwrite
from scipy.sparse import coo_matrix
import argparse
import os
import scanpy as sc

def fetch_umi(mtx_dir):
    adata = sc.read_10x_mtx(mtx_dir, var_names='gene_symbols', cache = False)
    return adata.to_df()

def main():
    parser = argparse.ArgumentParser(description='')
    parser.add_argument(f'--mtx_dir',help='', required=True)
    parser.add_argument(f'--sample_name',help='', required=True)
    parser.add_argument(f'--drop_bcs',help='', required=True)
    args = parser.parse_args()
    
    out_dir = f"./out_new_mtx/{sample_name}_10x_matrix/"
    genes = glob.glob(f"{args.mtx_dir}/genes.tsv")[0]
    os.system(f"mkdir -p {out_dir}")

    ft = fetch_umi(args.mtx_dir)

    with open(args.drop_bcs) as fd:
        drops = [i.strip("\n") for i in fd.readlines()]

    umi["barcode"] = list(umi.index)

    umi = umi[~umi["barcode"].isin(drops)]

    gene_id_dict = {}

    with open(genes) as fd:
        for record in fd.readlines():
            gene_id, gene_sy = record.split("\t")
            gene_id_dict.update({gene_sy.strip("\n"):gene_id.strip("\n")})
    # to mtx
    # write barcode
    barcodes = umi['barcode'].to_list()
    with open(f"{out_dir}/barcodes.tsv", "w") as fd:
        for bc in barcodes:
            fd.write(f"{bc}\n")

    # write gene
    target_gene_umi_count_df = umi.drop('barcode', axis = 1)
    genes = target_gene_umi_count_df.columns

    with open(f"{out_dir}/genes.tsv", "w") as fd:
        for ge in genes:
            fd.write(f"{gene_id_dict[ge]}\t{ge}\n")

    # write mtx
    matrix = target_gene_umi_count_df.T
    mtx = coo_matrix(matrix, matrix.shape, int)
    mmwrite(f'{out_dir}/matrix.mtx', mtx)

if __name__ == '__main__':
    main()
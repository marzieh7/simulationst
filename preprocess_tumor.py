import scanpy as sc
import os

def preprocess_tumor(input_file, output_file):
    adata = sc.read_h5ad(input_file)
    print("Original shape:", adata.shape)

    # 2. Quality Control (keep)
    adata.var['mt'] = adata.var_names.str.startswith('MT')
    sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)
    # sc.pl.violin(adata, ['n_genes_by_counts', 'total_counts', 'pct_counts_mt'], jitter=0.4, multi_panel=True)

    # data ident 22 , 23 remove for benchmark data set  , 18 , 25 , 20 for model
    adata = adata[ (adata.obs['orig.ident'] == 'Pre_P022_t')| 
              (adata.obs['orig.ident'] == 'Pre_P023_t'),:]

    adata = adata[adata.obs.n_genes_by_counts > 200, :]
    adata = adata[adata.obs.n_genes_by_counts < 7000, :]
    adata = adata[adata.obs.pct_counts_mt < 5, :]
    print("After filtering:", adata.shape)

    sc.pp.neighbors(adata)

    sc.tl.umap(adata)

    print(adata.obs['condition'])
    

    sc.pl.umap(adata, color="orig.ident", save="orig.iden_tumorsc.png")
    sc.pl.umap(adata, color="condition", save="condition_tumorsc.png")

    os.makedirs(os.path.dirname(output_file), exist_ok=True)
    adata.write(output_file)
    print(f"Saved preprocessed file at: {output_file}")
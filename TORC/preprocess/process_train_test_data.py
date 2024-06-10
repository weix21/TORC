import os, sys
import numpy as np
import pandas as pd
import scipy
import scanpy as sc
import anndata
from scipy.spatial import distance

## plot
import matplotlib.pyplot as plt

## ---- some functions for processing data
def process_adata(adata, min_genes=10, min_cells=10, celltype_label="cell.type"):
    '''Procedures for filtering single-cell data
       1. Filter low-quality cells and genes;
       2. Filter nonsense genes;
       3. Normalize and log-transform the data;
       4. Change all gene names into UPPER;
       5. Remove cells with no labels;
    '''
    adata.var_names=[i.upper() for i in list(adata.var_names)]#avod some genes having lower letter

    ## make names unique
    adata.var_names_make_unique()
    adata.obs_names_make_unique()

    #3 prefilter_specialgene: MT and ERCC  -> from ItClust package
    Gene1Pattern="ERCC"
    Gene2Pattern="MT-"
    id_tmp1=np.asarray([not str(name).startswith(Gene1Pattern) for name in adata.var_names],dtype=bool)
    id_tmp2=np.asarray([not str(name).startswith(Gene2Pattern) for name in adata.var_names],dtype=bool)
    id_tmp=np.logical_and(id_tmp1,id_tmp2)
    adata._inplace_subset_var(id_tmp)

    ## handel exception when there are not enough cells or genes after filtering
    if adata.shape[0] < 3 or adata.shape[1] < 3:
        return None

    #4 normalization,var.genes,log1p
    sc.pp.normalize_per_cell(adata,counts_per_cell_after=10000)  ## total count equal to the median of the counts_per_cell
    sc.pp.log1p(adata)

    ## cells with celltypes
    cells = adata.obs.dropna(subset=[celltype_label]).index.tolist()
    adata = adata[cells]
    return adata

def feature_selection_train_test(train_adata, test_adata, result_dir,
        gene_no=1000, select_on="test", select_method="F-test",
        min_genes=10, min_cells=10, celltype_label="cell.type"):
    '''Perform feature selection on train_adata and test_adata

    @ result_dir: when method is FEAST, for storing temp counts and FEAST features
    @ gene_no: top x numbers of features
    @ select_on: whether perform on test or train
    @ select_method: Seurat/FEAST/F-test
        FEAST(unsupervised)/F-test(supervised) based on FEAST implementation and
        can be only applied to training datasets because we do not know labels
        for test
    @ min_genes/min_cells: when processing the data, the minimun requirements
    '''

    ## if feature already exists
    feature_file = result_dir+os.sep+"features.txt"
    if os.path.exists(feature_file):
        ## filter cells/genes, etc
        train_adata = process_adata(train_adata, min_genes, min_cells)
        test_adata = process_adata(test_adata, min_genes, min_cells)

        with open(feature_file) as f:
            features = f.read().splitlines()
        print("Number of features:", len(features))

        features = set(train_adata.var_names.tolist()).intersection(set(features))
        features = set(test_adata.var_names.tolist()).intersection(set(features))
        features = list(features)
        features.sort()  ## for reproducibility

        ## order common genes in anndata
        train_adata = train_adata[:, features]
        test_adata = test_adata[:, features]
        return train_adata, test_adata

    ## filter cells/genes, etc
    train_adata = process_adata(train_adata, min_genes, min_cells)
    test_adata = process_adata(test_adata, min_genes, min_cells)

    ## handle None exception
    if train_adata is None or test_adata is None:
        return None, None

    ## select top 1000 HVGs from test
    if select_method == "Seurat":
        if select_on == "test":
            sc.pp.highly_variable_genes(test_adata, n_top_genes=gene_no, subset=True)
        if select_on == "train":
            sc.pp.highly_variable_genes(train_adata, n_top_genes=gene_no, subset=True)

    if "FEAST" == select_method or "F-test" == select_method:
        ## read features selected by FEAST
        if scipy.sparse.issparse(train_adata.X) or \
                isinstance(train_adata.X, pd.DataFrame):
            tmp_data = train_adata.X.toarray()
        else:
            tmp_data = train_adata.X

        ## calculate F-test
        cell_annots = train_adata.obs[celltype_label].tolist()
        uniq_celltypes = set(cell_annots)
        array_list = []
        for celltype in uniq_celltypes:
            idx = np.where(np.array(cell_annots) == celltype)[0].tolist()
            array_list.append(tmp_data[idx, :])
        F, p = scipy.stats.f_oneway(*array_list)
        F_updated = np.nan_to_num(F)
        sorted_idx = np.argsort(F_updated)[-gene_no:]
        feast_features = train_adata.var_names[sorted_idx].tolist()
        feast_features.sort()

        if select_on == "test":
            feast_genes = set(feast_features).intersection(set(test_adata.var_names.tolist()))
            test_adata = test_adata[:, list(feast_genes)]
        if select_on == "train":
            feast_genes = set(feast_features).intersection(set(train_adata.var_names.tolist()))
            train_adata = train_adata[:, list(feast_genes)]

    features = set(train_adata.var_names.tolist()).intersection(set(test_adata.var_names.tolist()))
    features = list(features)
    features.sort()  ## for reproducibility
    print("Number of features:", len(features))

    ## write features into file
    with open(result_dir+os.sep+"features.txt", 'w') as f:
        for feature in features:
            f.write("%s\n" % feature)

    ## write the numbers of features into file
    with open(result_dir+os.sep+"n_features.txt", 'w') as f:
        f.write("%s\n" % len(features))

    ## order common genes in anndata
    train_adata = train_adata[:, features]
    test_adata = test_adata[:, features]
    return train_adata, test_adata

def scale_and_visualize(train_adata, test_adata, result_dir, dr_seed=0, scale=True,
        plot=True, plot_elements=['dataset_batch', 'cell.type']):
    '''Scale data set and plot a dimension reduction on certain elements
    @dr_seed: seed for dimention reduction
    @plot_elements: dataset_batch, cell.type, or ind
    '''
    if scale:
        sc.pp.scale(train_adata, zero_center=True, max_value=6)
        X = test_adata.X.toarray()
        X -= train_adata.var["mean"].values
        X /= train_adata.var["std"].values
        X[X > 6] = 6
        test_adata.X = X
        del(X)

    ## plot train_adata and test_adata
    adata = train_adata.concatenate(test_adata,join='inner',
            batch_key="dataset_batch",batch_categories=["train","test"]) #inner join

    if plot:
        if "ind" in adata.obs.columns:
            adata.obs["ind"] = adata.obs["ind"].astype("category")
        plot_adata(adata, plot_elements, result_dir, dr_seed=dr_seed)

    ## set adata information to train_adata, test_adata
    train_adata = adata[adata.obs[adata.obs["dataset_batch"] == "train"].index.tolist()]
    test_adata = adata[adata.obs[adata.obs["dataset_batch"] == "test"].index.tolist()]
    return train_adata, test_adata

def load_adata(result_dir):
    '''If data already exists, load from disk
    '''
    if (os.path.exists(result_dir+os.sep+"train_adata.h5ad") and
            os.path.exists(result_dir+os.sep+"test_adata.h5ad")):
        train_adata = anndata.read_h5ad(result_dir+os.sep+"train_adata.h5ad")
        test_adata = anndata.read_h5ad(result_dir+os.sep+"test_adata.h5ad")
        return True, train_adata, test_adata
    else:
        return False, None, None

def save_adata(train_adata, test_adata, result_dir, write=True):
    '''Save data to disk
    '''
    if write:
        train_adata.write(result_dir+os.sep+"train_adata.h5ad")
        test_adata.write(result_dir+os.sep+"test_adata.h5ad")

def plot_adata(adata, columns, result_dir, dr_seed=0, prefix=""):
    '''Dimension reduction on adata and plot out features of interest

    @ adata: combined anndata of train and test
    @ columns: columns of interest from adata.obs
    @ result_dir: where to store the plots
    @ dr_seed: dimension reduction seed
    '''
    # do PCA first
    sc.tl.pca(adata, random_state=dr_seed)
    sc.tl.tsne(adata, use_rep="X_pca",
            learning_rate=300, perplexity=30, n_jobs=4, random_state=dr_seed)
    sc.pl.tsne(adata, color=columns)
    plt.savefig(result_dir+os.sep+prefix+"tSNE_cluster.png")

    ## do UMAP
    sc.pp.neighbors(adata, n_neighbors=20, use_rep="X_pca", random_state=dr_seed)
    sc.tl.umap(adata, random_state=dr_seed)
    sc.pl.umap(adata, color=columns)
    plt.savefig(result_dir+os.sep+prefix+"umap_cluster.png")

def process_pipeline(train_adata, test_adata, result_dir,
        gene_no=1000, select_on="test", select_method="Seurat",
        min_genes=10, min_cells=10):
    ''' A process pipeline integrating feature selection, center scaled the data
    '''
    ## feature selection
    train_adata, test_adata = feature_selection_train_test(train_adata, test_adata,
            result_dir, gene_no, select_on, select_method, min_genes, min_cells)

    ## if after feature selection, one of them is None, then train and test to None
    if train_adata is None or test_adata is None:
        return None, None

    ## scale and analyze
    train_adata, test_adata = scale_and_visualize(train_adata, test_adata, result_dir,
            plot=False)
    return train_adata, test_adata

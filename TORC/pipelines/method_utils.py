import os
import math
import sys
import logging
import anndata
import numpy as np
import pandas as pd
import scanpy as sc
import scipy
import random
import tensorflow as tf


import matplotlib.pyplot as plt
plt.rcParams.update({'font.size': 18})

from sklearn.preprocessing import OneHotEncoder
from sklearn import metrics

from typing import TypeVar
A = TypeVar('anndata')  ## generic for anndata
ENC = TypeVar('OneHotEncoder')

from models.MLP import MLP
from preprocess.process_data_train_test import _process_initial
from preprocess.process_data_train_test import _process_reconstruct

logger = logging.getLogger(__name__)

# RANDOM_SEED = 1993
# random.seed(RANDOM_SEED)
# np.random.seed(RANDOM_SEED)

MLP_DIMS = [64, 16]
Teacher_DIMS = [64, 16]
Student_DIMS = [64, 16]
BATCH_SIZE = 32
Celltype_COLUMN = "cell.type"
PredCelltype_COLUMN = "pred_celltype"
ENTROPY_QUANTILE = 0.4  ## how many target cells could be used for expanding reference


### === load anndata
def _load_adata(args=None):
    train_adata, test_adata = None, None
    ## === TORC 
    if args.torc_step == "TORC_initial":
        train_adata, test_adata = \
            _process_initial(args)
    if args.torc_step == "TORC_reconstruct":
        train_adata, test_adata = \
            _process_reconstruct(args)

    return train_adata, test_adata

def _process_loaded_data(train_adata, test_adata, result_dir,
        args=None, save_raw=False):
     
    if args is None:
        sys.exit("Error: Please check your argument parser object!")

    if save_raw:
        train_adata.layers["counts"] = train_adata.X.copy()
        test_adata.layers["counts"] = test_adata.X.copy()


    ## Process data
    train_adata = _process_adata(train_adata, process_type='train')
    test_adata = _process_adata(test_adata, process_type='test')

    ## feature selection
    train_adata = _select_feature(train_adata, fs_method = args.select_method, num_features = args.n_features)
    features = set(train_adata.var_names.tolist()).intersection(set(test_adata.var_names.tolist()))
    features = list(features)
    train_adata = train_adata[:, features]
    test_adata = test_adata[:, features]

    # scale data
    if args.scale:
        train_adata = _scale_data(train_adata)
        test_adata = _scale_data(test_adata)

    return train_adata, test_adata



def _process_adata(adata, process_type='train', celltype_label='cell.type'):
    '''Procedures for filtering single-cell gene scale data (can be gene expression, or gene scores)
       1. Filter nonsense genes;
       2. Normalize and log-transform the data;
       3. Remove cells with no labels;
    '''
    adata.var_names=[i.upper() for i in list(adata.var_names)] #avoid some genes having lower letter

    ## make names unique after removing
    adata.var_names_make_unique()
    adata.obs_names_make_unique()

    #prefilter_specialgene: MT and ERCC  -> refered from ItClust package
    Gene1Pattern="ERCC"
    Gene2Pattern="MT-"
    id_tmp1=np.asarray([not str(name).startswith(Gene1Pattern) for name in adata.var_names],dtype=bool)
    id_tmp2=np.asarray([not str(name).startswith(Gene2Pattern) for name in adata.var_names],dtype=bool)
    id_tmp=np.logical_and(id_tmp1,id_tmp2)
    adata._inplace_subset_var(id_tmp)

    ## handel exception when there are not enough cells or genes after filtering
    if adata.shape[0] < 3 or adata.shape[1] < 3:
        sys.exit("Error: too few genes or cells left to continue..")

    ## normalization,var.genes,log1p
    sc.pp.normalize_per_cell(adata, counts_per_cell_after=10000, min_counts=0)
    sc.pp.log1p(adata)

    ## cells with celltypes
    if process_type == 'train':
        cells = adata.obs.dropna(subset=[celltype_label]).index.tolist()
        adata = adata[cells]
    return adata

def _select_feature(adata: A, fs_method = "F-test", num_features: int = 1000) -> A:
    '''Select features
    ---
    Input:
        - anndata
        - fs_method: F-test / noFS / Seurat
    '''
    ## Feature selection
    if fs_method == "noFS":
        logger.info("Feature selection will not be performed.\n")
        return adata

    if fs_method == "F-test":
        print("Use F-test to select features.\n")
        if scipy.sparse.issparse(adata.X) or \
                isinstance(adata.X, pd.DataFrame):
            tmp_data = adata.X.toarray()
        else:
            tmp_data = adata.X

        ## calculate F-test
        cell_annots = adata.obs[Celltype_COLUMN].tolist()
        uniq_celltypes = set(cell_annots)
        array_list = []
        for celltype in uniq_celltypes:
            idx = np.where(np.array(cell_annots) == celltype)[0].tolist()
            array_list.append(tmp_data[idx, :])
        F, p = scipy.stats.f_oneway(*array_list)
        F_updated = np.nan_to_num(F)
        sorted_idx = np.argsort(F_updated)[-num_features:]
        features = adata.var_names[sorted_idx].tolist()
        features.sort()
        adata = adata[:, features]

    if fs_method == "Seurat":
        print("Use Seurat in scanpy to select features.\n")
        sc.pp.highly_variable_genes(adata, n_top_genes=num_features, subset=True)

    return adata


def _scale_data(adata):
    '''Center scale
    '''
    sc.pp.scale(adata, zero_center=True, max_value=6)
    return adata

def _visualize_data(adata, output_dir, color_columns=["celltype"],
        reduction="tSNE", prefix="data"):
    '''Visualize data
    ---
    Input:
        - reduction: tSNE or UMAP
        - color_columns: plot on categories
    '''
    sc.tl.pca(adata, random_state=RANDOM_SEED)

    if reduction == "tSNE":
        sc.tl.tsne(adata, use_rep="X_pca",
            learning_rate=300, perplexity=30, n_jobs=1, random_state=RANDOM_SEED)
        sc.pl.tsne(adata, color=color_columns)
        plt.tight_layout()
        plt.savefig(output_dir+os.sep+prefix+"tSNE_cluster.png")
    if reduction == "UMAP":
        sc.pp.neighbors(adata, n_neighbors=20, use_rep="X_pca", random_state=RANDOM_SEED)
        sc.tl.umap(adata, random_state=RANDOM_SEED)
        sc.pl.umap(adata, color=color_columns)
        plt.tight_layout()
        plt.savefig(output_dir+os.sep+prefix+"umap_cluster.png")

def _save_adata(adata, output_dir, prefix=""):
    '''Save anndata as h5ad
    '''
    adata.write(output_dir+os.sep+prefix+'adata.h5ad')


def _prob_to_label(y_pred: np.ndarray, encoders: dict) -> list:
    '''Turn predicted probabilites to labels
    ---
    Input:
        - y_pred: Predicted probabilities
        - encoders: dictionary with mapping information
    ---
    Output:
        - a list containing predicted cell types
    '''
    pred_labels = y_pred.argmax(1)
    pred_celltypes = [encoders[label] for label in pred_labels]
    print("=== Predicted celltypes: ", set(pred_celltypes))
    return pred_celltypes

def _label_to_onehot(labels: list, encoders:dict) -> np.ndarray:
    '''Turn predicted labels to onehot encoder
    ---
    Input:
        - labels: the input predicted cell types
        - encoders: dictionary with mapping information
    '''
    inv_enc = {v: k for k, v in encoders.items()}
    onehot_arr = np.zeros((len(labels), len(encoders)))
    pred_idx = [inv_enc[l] for l in labels]
    onehot_arr[np.arange(len(labels)), pred_idx] = 1
    return onehot_arr


def _extract_adata(adata: A) -> np.ndarray:
    '''Extract adata.X to a numpy array
    ---
    Output:
         - matrix in np.ndarray format
    '''
    if scipy.sparse.issparse(adata.X) or isinstance(adata.X, pd.DataFrame) or isinstance(adata.X, anndata._core.views.ArrayView):
        X = adata.X.toarray()
    else:
        X = adata.X
    return X

def _init_MLP(x_train, y_train, dims=[64, 16], seed=0):
    '''Initialize MLP model based on input data
    '''
    mlp = MLP(dims)
    mlp.input_shape = (x_train.shape[1], )
    #mlp.n_classes = len(set(y_train.argmax(1)))
    mlp.n_classes = y_train.shape[1]
    mlp.random_state = seed
    mlp.init_MLP_model()  ## init the model
    return mlp

def _select_confident_cells(adata, celltype_col):
    '''Select low entropy cells from each predicted cell type
    ---
    Input:
        - adata: anndata object
        - celltype_col: the column indicator
    '''
    low_entropy_cells = []
    for celltype in set(adata.obs[celltype_col]):
        celltype_df = adata.obs[adata.obs[celltype_col] == celltype]
        entropy_cutoff = np.quantile(celltype_df['entropy'], q=ENTROPY_QUANTILE)
        ## change to < instead of <= to deal with ties
        cells = celltype_df.index[np.where(celltype_df['entropy'] <= entropy_cutoff)[0]].tolist()
        num_cells = math.ceil(ENTROPY_QUANTILE*celltype_df.shape[0])
        if len(cells) > num_cells:
            random.seed(RANDOM_SEED)
            selected_cells = random.sample(cells, num_cells)
        else:
            selected_cells = cells
        low_entropy_cells.extend(selected_cells)
    high_entropy_cells = list(set(adata.obs_names) - set(low_entropy_cells))
    adata.obs.loc[low_entropy_cells, 'entropy_status'] = "low"
    adata.obs.loc[high_entropy_cells, 'entropy_status'] = "high"
    return adata


def _run_pipeline(train_adata, test_adata, result_dir, rseed=1234):

    RANDOM_SEED = rseed
    random.seed(RANDOM_SEED)
    np.random.seed(RANDOM_SEED)

    ## Train model
    x_train = _extract_adata(train_adata)
    enc = OneHotEncoder(handle_unknown='ignore')
    y_train = enc.fit_transform(train_adata.obs[["cell.type"]]).toarray()
    logger.debug("Categories information: ", enc.categories_[0])

    mlp = _init_MLP(x_train, y_train, dims=MLP_DIMS,
            seed=RANDOM_SEED)
    mlp.compile()
    mlp.fit(x_train, y_train)

    with open(result_dir+os.sep+"onehot_encoder.txt", 'w') as f:
        for idx, cat in enumerate(enc.categories_[0]):
            f.write('%d:%s\n' % (idx, cat))

    model = mlp.model
    encoders = {}
    with open(result_dir+os.sep+"onehot_encoder.txt") as f:
        for line in f:
            line_info = line.strip().split(':')
            encoders[int(line_info[0])] = line_info[1]


    ## Naive MLP
    test_data_mat = _extract_adata(test_adata)
    y_pred = tf.nn.softmax(model.predict(test_data_mat)).numpy()
    pred_celltypes = _prob_to_label(y_pred, encoders)
    test_adata.obs[PredCelltype_COLUMN] = pred_celltypes

    y_avg_pred = np.average(y_pred,axis=0)
    df = pd.DataFrame(y_avg_pred)
    df.columns = ["prob"]
    df.index = enc.categories_[0]
    df.to_csv(result_dir+os.sep+'avg_post_prob.csv')

    entropy = [-np.nansum(y_pred[i]*np.log(y_pred[i])) for i in range(y_pred.shape[0])]
    test_adata.obs['entropy'] = entropy
    test_adata = _select_confident_cells(
                test_adata, celltype_col=PredCelltype_COLUMN)

    ## select certain columns and store to the file
    test_adata.obs.to_csv(result_dir+os.sep+"predicted_obs.csv")

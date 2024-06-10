'''
Load data for running celltyping experiments
'''
import os, sys
from preprocess.process_train_test_data import *

from preprocess.process_PBMC_train_test import *

#Distance_RSCRIPT_PATH = "Distance.R"

### === process loaded data
def process_loaded_data(train_adata, test_adata, result_dir,
        args=None, scale=True, plot=False,
        save_raw=False, save_data=False):

    if args is None:
        sys.exit("Error: Please check your argument parser object!")

    #curate for common cell types
    common_celltypes = args.common_ct
    common_celltypes = common_celltypes.split(',')
    train_cells = train_adata.obs.loc[train_adata.obs["cell.type"].isin(common_celltypes)].index
    test_cells = test_adata.obs.loc[test_adata.obs["cell.type"].isin(common_celltypes)].index
    train_adata = train_adata[train_cells]
    test_adata = test_adata[test_cells]
    print("train_adata: \n", set(train_adata.obs["cell.type"]))
    print("test_adata: \n", set(test_adata.obs["cell.type"]))


    if save_raw:
        train_adata.layers["counts"] = train_adata.X.copy()
        test_adata.layers["counts"] = test_adata.X.copy()

    ## feature selection
    train_adata, test_adata = feature_selection_train_test(train_adata, test_adata,
            result_dir, args.n_features, args.select_on, args.select_method)

    if save_data:
        save_adata(train_adata, test_adata, result_dir)

    # scale and analze, no need to scale for calculating distance
    if args.scale:
        print("Do scale")
        train_adata, test_adata = scale_and_visualize(train_adata, test_adata,
            result_dir, scale=scale, plot=plot)

    return train_adata, test_adata

### === load PBMC anndata
def load_PBMC_adata(result_dir, args=None, scale=True, plot=False):
    train_adata, test_adata = None, None
    ## === PBMC datasets
    if args.data_source == "PBMC_all":
        train_adata, test_adata = \
            process_all(result_dir,args)
    if args.data_source == "PBMC_sample":
        train_adata, test_adata = \
            process_sample(result_dir,args)
    return train_adata, test_adata

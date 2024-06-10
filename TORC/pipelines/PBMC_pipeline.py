'''
Configuration generation for running PBMC datasets
'''

import os, argparse
import io
import time

from pipelines import method_utils, dataloading_utils
from preprocess.process_train_test_data import *

sys.stdout = io.TextIOWrapper(sys.stdout.buffer,encoding='utf-8')

if __name__ == "__main__":

    ## parse arguments
    parser = argparse.ArgumentParser(description="Celltyping pipeline.")
    parser.add_argument('data_source', help="Load which dataset",
        choices=["PBMC_all","PBMC_sample"])

    parser.add_argument('-m', '--method', help="Run which method",
        choices=['MLP', 'SVM_RBF', 'SVM_linear'],
        required=True)
    parser.add_argument('--select_on', help="Feature selection on train or test, or None of them",
        choices=['train', 'test'])
    parser.add_argument('--select_method', help="Feature selection method, Seurat/FEAST or None",
            choices=['Seurat', 'FEAST', 'F-test'])
    parser.add_argument('--n_features', help="Number of features selected",
            default=1000, type=int)
    parser.add_argument('--train', help="Specify which as train")
    parser.add_argument('--test', help="Specify which as test")
    parser.add_argument('--sample_seed', help="Downsample seed in combined individual effect",
            default=0, type=int)
    parser.add_argument('--distance', help="Whether to calculate the distance",
            default=False)
    parser.add_argument('--pip_dir', help="Result Dir",
            required=True)
    parser.add_argument('--NMLP_dir', help="Naive MLP Result Dir")
    parser.add_argument('--data_dir', help="Target data Dir")
    parser.add_argument('--ref_dir', help="Ref data Dir")
    parser.add_argument('--sample_method', help="Define Sampling method")
    parser.add_argument('--scale', help="Define if scale", default=False)
    parser.add_argument('--file_name', help="Define file name")
    parser.add_argument('--rseed', help="Define the random seed used", default=1234, type=int)
    parser.add_argument('--common_ct',help="Define the celltypes for analyses",required=True)
    parser.add_argument('--sample_size',help="Roughly define the sample size for analysis", default=1000, type=int)
    parser.add_argument('--expand_ref', help="Define if expanding the reference")
    parser.add_argument('--Entropy_dir', help="Entropy Status Dir")

    args = parser.parse_args()

    start = time.time()

    result_dir = args.pip_dir
    os.makedirs(result_dir, exist_ok=True)

    load_ind, train_adata, test_adata = load_adata(result_dir)
    if not load_ind:
        train_adata, test_adata = dataloading_utils.load_PBMC_adata(
            result_dir, args=args)

        train_adata, test_adata = dataloading_utils.process_loaded_data(
                train_adata, test_adata, result_dir, args=args)
        print("Train anndata: \n", train_adata)
        print("Test anndata: \n", test_adata)

    test_adata = method_utils._run_pipeline(train_adata, test_adata, result_dir, rseed=args.rseed)

    end = time.time()

    runTime = end - start

    with open(result_dir+os.sep+"Runtime.txt", 'w') as f:
        f.write("Time:%s\n" % str(runTime))

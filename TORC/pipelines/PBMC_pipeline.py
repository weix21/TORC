import os, argparse

from pipelines import method_utils

if __name__ == "__main__":

    ## parse arguments
    parser = argparse.ArgumentParser(description="Celltyping pipeline.")
    parser.add_argument('--torc_step', help="Which step of TORC",
        choices=["TORC_initial","TORC_reconstruct"])
    parser.add_argument('--select_method', help="Feature selection method, Seurat, F-test or None",
            choices=['noFS', 'F-test', 'Seurat'])
    parser.add_argument('--n_features', help="Number of features selected",
            default=1000, type=int)
    parser.add_argument('--res_dir', help="Result Dir",
            required=True)
    parser.add_argument('--NMLP_dir', help="Initial MLP Result Dir")
    parser.add_argument('--data_dir', help="Target data Dir")
    parser.add_argument('--ref_dir', help="Ref data Dir")
    parser.add_argument('--train_sample', help="Train sample ID")
    parser.add_argument('--test_sample', help="Test sample ID")
    parser.add_argument('--sample_method', help="Define Sampling method",
                        choices=['Prop', 'Prob'])
    parser.add_argument('--rseed', help="Define the random seed used", default=1234, type=int)
    parser.add_argument('--sample_size',help="Roughly define the sample size for analysis", default=1000, type=int)
    parser.add_argument('--expand_ref', help="Define if expanding the reference")
    parser.add_argument('--Entropy_dir', help="Entropy Status Dir")

    args = parser.parse_args()

    result_dir = args.res_dir
    os.makedirs(result_dir, exist_ok=True)
    
    train_adata, test_adata = method_utils._load_adata(
            args=args)

    train_adata, test_adata = method_utils._process_loaded_data(
                train_adata, test_adata, result_dir, args=args)
    
    print("Train anndata: \n", train_adata)
    print("Test anndata: \n", test_adata)

    test_adata = method_utils._run_pipeline(train_adata, test_adata, result_dir, rseed=args.rseed)


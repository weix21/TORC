import anndata
import random
import copy
import pandas as pd
import math
import numpy as np

from preprocess import load_PBMC_data


def process_all(result_dir, args,
        celltype_gran=1):
    ''' Process individual data of PBMC batch2
    @data_dir: where PBMC batch2 data stroes
    @result_dir: where to store PCA/tSNE/UMAP result
    @input1/input2: can be batch1_indID/batch2_indID
    @celltype_gran: granularity of cell types, 0: major cell types, 1:sub-celltypes
        According to the dataset, the give cell types are sub-cell types
    '''
    ## split input1
    if args.train is None:
        train_adata = anndata.read_h5ad(args.ref_dir)
    else:
        print("Sepcify samples are used")
        adata  = anndata.read_h5ad(args.ref_dir)
        ind_cells = adata.obs[adata.obs["ind"].isin(args.train.split(','))].index
        train_adata = adata[ind_cells]

    if args.test is None:
        test_adata  = anndata.read_h5ad(args.data_dir)
    else:
        print("Sepcify samples are used")
        adata  = anndata.read_h5ad(args.data_dir)
        ind_cells = adata.obs[adata.obs["ind"].isin(args.test.split(','))].index
        test_adata = adata[ind_cells]
        del(adata)

    print("Train anndata: \n", train_adata)
    print("Test anndata: \n", test_adata)

    ## curate given sub-cell types to major cell types
    train_adata = load_PBMC_data.curate_PBMC_demulx_celltypes(train_adata, celltype_gran)
    print("train_adata: \n", set(train_adata.obs["cell.type"]))
    test_adata = load_PBMC_data.curate_PBMC_demulx_celltypes(test_adata, celltype_gran)
    print("test_adata: \n", set(test_adata.obs["cell.type"]))
    return train_adata, test_adata




def process_sample(result_dir, args,
        celltype_gran=1):
    ''' Process individual data of PBMC batch2
    @data_dir: where PBMC batch2 data stroes
    @result_dir: where to store PCA/tSNE/UMAP result
    @input1/input2: can be batch1_indID/batch2_indID
    @celltype_gran: granularity of cell types, 0: major cell types, 1:sub-celltypes
        According to the dataset, the give cell types are sub-cell types
    '''

    RANDOM_SEED = args.rseed
    random.seed(RANDOM_SEED)
    sample_cells = []
    metadata = pd.read_csv(args.NMLP_dir,index_col=0)

    if 'scNym' in metadata.columns:
        metadata['firstround_pred_celltype'] = metadata['scNym']
    elif 'scanvi' in metadata.columns:
        metadata['firstround_pred_celltype'] = metadata['scanvi']
    elif 'predictions' in metadata.columns:
        metadata['firstround_pred_celltype'] = metadata['predictions']

    if args.train is None:
        all_adata = anndata.read_h5ad(args.ref_dir)
    else:
        adata  = anndata.read_h5ad(args.ref_dir)
        ind_cells = adata.obs[adata.obs["ind"].isin(args.train.split(','))].index
        all_adata = adata[ind_cells]

    if args.test is None:
        test_adata  = anndata.read_h5ad(args.data_dir)
    else:
        adata  = anndata.read_h5ad(args.data_dir)
        ind_cells = adata.obs[adata.obs["ind"].isin(args.test.split(','))].index
        test_adata = adata[ind_cells]
        del(adata)

    test_adata = load_PBMC_data.curate_PBMC_demulx_celltypes(test_adata, celltype_gran=1)
    all_adata = load_PBMC_data.curate_PBMC_demulx_celltypes(all_adata, celltype_gran=1)

    if (args.expand_ref == "Expand_v1"):
        print("Expand reference")
        entropy_data = pd.read_csv(args.Entropy_dir,index_col=0)
        if 'scNym' in entropy_data.columns:
            entropy_data['firstround_pred_celltype'] = entropy_data['scNym']
        elif 'scanvi' in entropy_data.columns:
            entropy_data['firstround_pred_celltype'] = entropy_data['scanvi']
        elif 'predictions' in entropy_data.columns:
            entropy_data['firstround_pred_celltype'] = entropy_data['predictions']
        low_entropy_cells = entropy_data.loc[entropy_data['entropy_status'] == 'low']["barcode"].tolist()
        high_entropy_cells = entropy_data.loc[entropy_data['entropy_status'] == 'high']["barcode"].tolist()
        test_ref_adata = test_adata[low_entropy_cells]
        test_ref_adata.obs["cell.type"] = entropy_data.loc[entropy_data['entropy_status'] == 'low']["firstround_pred_celltype"].tolist()

        if test_ref_adata.n_obs > all_adata.n_obs:
            test_ref_adata_number = all_adata.n_obs
            test_ref_proptionDict = dict(pd.Series(test_ref_adata.obs["cell.type"]).value_counts(normalize=True))
            test_ref_sample_cells = []
            for celltype in sorted(set(test_ref_proptionDict)):
                cells=[]
                cell_n = math.ceil(test_ref_proptionDict[celltype]*test_ref_adata_number)
                celltype_df = test_ref_adata.obs[test_ref_adata.obs["cell.type"] == celltype]
                if cell_n < 50:
                    cell_n = 50
                print(celltype)
                print(cell_n)
                random.seed(RANDOM_SEED)
                if celltype_df.shape[0] >= cell_n:
                    cells.extend(random.sample(celltype_df.index.tolist(), cell_n))
                    print("n_cells \n", len(cells))
                else:
                    cells.extend(celltype_df.index.tolist())
                    print("n_cells \n", len(cells))
                    cells.extend(random.choices(celltype_df.index.tolist(), k=cell_n-len(cells)))
                    print("n_cells \n", len(cells))
                test_ref_sample_cells.extend(cells)
            test_ref_adata = test_ref_adata[test_ref_sample_cells]
            test_ref_adata.obs_names_make_unique(join="-")
            test_ref_adata.obs["barcode"] = test_ref_adata.obs_names

        print("test_ref_adata_size: \n", test_ref_adata.n_obs)
        all_adata = all_adata.concatenate(test_ref_adata,join='inner',
                batch_key="dataset_batch",batch_categories=["train","test_ref"]) #inner join
        print(all_adata)

    common_celltypes = args.common_ct
    common_celltypes = common_celltypes.split(',')
    if(args.sample_method == "Prop"):
        print("Based on Celltype Proportion")
        if(set(pd.Series(metadata["firstround_pred_celltype"])) == set(common_celltypes)):
            proptionDict = dict(pd.Series(metadata["firstround_pred_celltype"]).value_counts(normalize=True))
        else:
            print("Missing celltypes")
            proptionDict = dict(pd.Series(metadata["firstround_pred_celltype"]).value_counts()+1)
            for i in common_celltypes:
                if(i not in proptionDict):
                    proptionDict[i]=1
            propsum = sum([proptionDict[i] for i in proptionDict])
            for x in proptionDict:
                proptionDict[x] = proptionDict[x]/propsum

    if(args.sample_method == "Prob"):
        print("Based on Posterior prob")
        if(set(pd.Series(metadata.index)) == set(common_celltypes)):
            proptionDict = dict(pd.Series(metadata["prob"]))
            propsum = sum([proptionDict[i] for i in proptionDict])
            for x in proptionDict:
                proptionDict[x] = proptionDict[x]/propsum
        else:
            print("Missing celltypes")
            proptionDict = dict(pd.Series(metadata["prob"])+0.0001)
            for i in set(common_celltypes):
                if(i not in proptionDict):
                    proptionDict[i]=0.0001
            propsum = sum([proptionDict[i] for i in proptionDict])
            for x in proptionDict:
                proptionDict[x] = proptionDict[x]/propsum

    if(args.sample_method == "GT"):
        print("Based on Ground Truth")
        if(set(pd.Series(metadata["cell.type"])) == set(common_celltypes)):
            proptionDict = dict(pd.Series(metadata["cell.type"]).value_counts(normalize=True))
        else:
            print("Missing celltypes")
            proptionDict = dict(pd.Series(metadata["cell.type"]).value_counts()+1)
            for i in common_celltypes:
                if(i not in proptionDict):
                    proptionDict[i]=1
            propsum = sum([proptionDict[i] for i in proptionDict])
            for x in proptionDict:
                proptionDict[x] = proptionDict[x]/propsum

    print(proptionDict)
    print(len(set(proptionDict)))
    sample_cells_number = args.sample_size*len(set(proptionDict))

    for celltype in sorted(set(proptionDict)):
        cells=[]
        cell_n = math.ceil(proptionDict[celltype]*sample_cells_number)
        celltype_df = all_adata.obs[all_adata.obs["cell.type"] == celltype]
        if cell_n < 50:
            cell_n = 50
        print(celltype)
        print(cell_n)
        random.seed(RANDOM_SEED)
        if celltype_df.shape[0] >= cell_n:
            cells.extend(random.sample(celltype_df.index.tolist(), cell_n))
            print("n_cells \n", len(cells))
        else:
            cells.extend(celltype_df.index.tolist())
            print("n_cells \n", len(cells))
            cells.extend(random.choices(celltype_df.index.tolist(), k=cell_n-len(cells)))
            print("n_cells \n", len(cells))
        sample_cells.extend(cells)

    train_adata = all_adata[sample_cells]
    train_adata.obs_names_make_unique(join="-")
    train_adata.obs["barcode"] = train_adata.obs_names


    print("train_adata: \n", set(train_adata.obs["cell.type"]))
    print("test_adata: \n", set(test_adata.obs["cell.type"]))
    print(proptionDict)
    print(dict(pd.Series(train_adata.obs["cell.type"]).value_counts(normalize=True)))
    return train_adata, test_adata

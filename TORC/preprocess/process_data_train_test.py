import anndata
import random
import copy
import pandas as pd
import math
import numpy as np


PredCelltype_COLUMN = "pred_celltype"
Celltype_COLUMN = "cell.type"


def _process_initial(args):
    ''' Load Reference data and target data
    '''
    if args.train_sample is None:
        train_adata = anndata.read_h5ad(args.ref_dir)
    else:
        print("Sepcify samples are used")
        adata  = anndata.read_h5ad(args.ref_dir)
        ind_cells = adata.obs[adata.obs["ind"].isin(args.train.split(','))].index
        train_adata = adata[ind_cells]

    if args.test_sample is None:
        test_adata  = anndata.read_h5ad(args.data_dir)
    else:
        print("Sepcify samples are used")
        adata  = anndata.read_h5ad(args.data_dir)
        ind_cells = adata.obs[adata.obs["ind"].isin(args.test.split(','))].index
        test_adata = adata[ind_cells]
        del(adata)

    print("Train anndata: \n", train_adata)
    print("Test anndata: \n", test_adata)

    return train_adata, test_adata




def _process_reconstruct(args):
    ''' Load Reference data and target data, and perform reference construction
    '''

    RANDOM_SEED = args.rseed
    random.seed(RANDOM_SEED)
    sample_cells = []
    metadata = pd.read_csv(args.NMLP_dir,index_col=0)

    if args.train_sample is None:
        all_adata = anndata.read_h5ad(args.ref_dir)
    else:
        adata  = anndata.read_h5ad(args.ref_dir)
        ind_cells = adata.obs[adata.obs["ind"].isin(args.train.split(','))].index
        all_adata = adata[ind_cells]

    if args.test_sample is None:
        test_adata  = anndata.read_h5ad(args.data_dir)
    else:
        adata  = anndata.read_h5ad(args.data_dir)
        ind_cells = adata.obs[adata.obs["ind"].isin(args.test.split(','))].index
        test_adata = adata[ind_cells]
        del(adata)

    if (args.expand_ref == "Expand_v1"):
        print("Expand reference")
        entropy_data = pd.read_csv(args.Entropy_dir,index_col=0)

        low_entropy_cells = entropy_data.loc[entropy_data['entropy_status'] == 'low']["barcode"].tolist()
        high_entropy_cells = entropy_data.loc[entropy_data['entropy_status'] == 'high']["barcode"].tolist()
        test_ref_adata = test_adata[low_entropy_cells]
        test_ref_adata.obs[Celltype_COLUMN] = entropy_data.loc[entropy_data['entropy_status'] == 'low'][PredCelltype_COLUMN].tolist()

        if test_ref_adata.n_obs > all_adata.n_obs:
            test_ref_adata_number = all_adata.n_obs
            test_ref_proptionDict = dict(pd.Series(test_ref_adata.obs[Celltype_COLUMN]).value_counts(normalize=True))
            test_ref_sample_cells = []
            for celltype in sorted(set(test_ref_proptionDict)):
                cells=[]
                cell_n = math.ceil(test_ref_proptionDict[celltype]*test_ref_adata_number)
                celltype_df = test_ref_adata.obs[test_ref_adata.obs[Celltype_COLUMN] == celltype]
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





    common_celltypes = set(all_adata.obs[Celltype_COLUMN])

    if(args.sample_method == "Prop"):
        print("Based on Celltype Proportion")
        if(set(pd.Series(metadata[PredCelltype_COLUMN])) == set(common_celltypes)):
            proptionDict = dict(pd.Series(metadata[PredCelltype_COLUMN]).value_counts(normalize=True))
        else:
            print("Missing celltypes")
            proptionDict = dict(pd.Series(metadata[PredCelltype_COLUMN]).value_counts()+1)
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

    return train_adata, test_adata


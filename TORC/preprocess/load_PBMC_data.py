import os
import anndata
import numpy as np
import pandas as pd

## --- curate PBMC demultiplex cell types to two granularity
def curate_PBMC_demulx_celltypes(adata, celltype_gran=1,
        celltype_label="cell.type"):
    '''Curate PBMC demulx cell types into two-starge prediction

    @adata: anndata object
    @celltype_gran: default 1 because the given cell.type are sub-cell types;
        if 0, then curate to major cell types
    @celltype_label: column indicator for cell type in adata.obs
    '''

    celltype_categories = {'B cells': ['B cells', "B", 'B_exhausted', 'B_immature','B_malignant','B_naive','B_non-switched_memory','B_switched_memory'],
            'CD4 T cells': ['CD4 T cells', 'CD4','CD4.CM','CD4.EM','CD4.IL22','CD4.Naive','CD4.Prolif','CD4.Tfh','CD4.Th1','CD4.Th17','CD4.Th2','Memory T cells','Naive T cells','Regulatory T cells','CD4+ Naive','CD4+ Effector Memory','CD4+ Memory','CD4+ Effector-GZMK','CD4+ Effector-GNLY'],
            'CD8 T cells': ['CD8', 'CD8 T cells','CD8.EM','CD8.Naive','CD8.Prolif','CD8.TE','Cytotoxic T cells','Cytotoxic T cell', 'Naive Cytotoxic T cells','CD8+ Effector-GZMK','CD8+ Effector-GNLY','CD8+ Naive'],
            'NK cells': ['NK cells', 'NK','NK_16hi','NK_56hi','NK_prolif','CD56 NK cells'],
            'Dendritic cells': ['Dendritic cells', 'DC','DC1','DC2','DC3','DC_prolif','pDC','ASDC','Plasmacytoid Dendritic cells ','Dendritic cells ','Mono DCs'],
            'Monocytes': ['CD14+ Monocytes', 'FCGR3A+ Monocytes', 'CD16+ monocyte cell', 'Mono', 'Monocytes','C1_CD16_mono','CD14_mono','CD16_mono','CD83_CD14_mono','Mono_prolif','CD14+ Monocytes', 'FCGR3A+ Monocytes', 'CD16+ monocyte cell','Monocytes','CD14+ Mono','CD16+ Mono','CD4+ Naive'],
            'Megakaryocytes': ['Megakaryocytes','Mega'],
            'Neu': ['Neu'],
            'Macro': ['Macro'],
            'Plasma': ['Plasma']
            }
    major_celltypes = []
    for celltype in adata.obs["cell.type"].tolist():
        flag = False
        for major_celltype, sub_celltypes in celltype_categories.items():
            if celltype in sub_celltypes:
                major_celltypes.append(major_celltype)

                ## if found
                flag = True
                break
        if False == flag:
            major_celltypes.append('unknown')
    ###Notice here
    adata.obs['cell.type'] = major_celltypes
    #adata.obs['majortypes'] = major_celltypes


    if 0 == celltype_gran:
        adata.obs.rename(columns={celltype_label: "subtypes",
                                    "majortypes": celltype_label},
                        inplace=True)
    return adata

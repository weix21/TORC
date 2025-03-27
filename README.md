# TORC:Target-Oriented Reference Construction for supervised cell-type identification in scRNA-seq

## Overview
TORC aims to construct a reference suitable for a given target data. It first uses an off-the-shelf supervised method to label the target cells, from which TORC first estimates cell-type composition in the target. TORC then add target cells with high-confidence labels to the reference to form an expanded reference pool. Then TORC resamples from the pool to construct a new reference according to the target composition. The reconstructed reference is used to build the final classifier. 

<p align="center">
  <img src="https://github.com/user-attachments/assets/bcb3fb94-a70e-41ca-862e-5ac210d605a6" />
</p>

Since TORC is a broadly applicable strategy for most supervised methods, this GitHub repository demonstrates how we apply TORC to our designed MLP as an example.

All the copyrights are explained by Xin Wei <xin_wei@brown.edu> from [Dr. Zhijin Wu's lab](https://vivo.brown.edu/display/zjwu) and [Dr. Hao Wu's lab](http://www.haowulab.org/). 
Any TORC-related questions should be posted to the GitHub Issue section of `TORC`
homepage at https://github.com/weix21/TORC/issues.

## Software requirements

The example requires the following:
- python=3.8
- numpy==1.19.2
- pandas==1.4.2
- scipy==1.8.1
- anndata==0.7.4
- scanpy==1.8.2
- matplotlib==3.5.2
- scikit-learn==1.1.1
- tensorflow==2.4.1
- keras==2.9.0
- statsmodels==0.13.2


## Predict cell types

### To get the initial predicted labels:
```
python -m pipelines.PBMC_pipeline --torc_step TORC_initial --select_method $fs_method --res_dir $RES_DIR --data_dir $DATA_DIR --ref_dir $REF_DIR --scale

```

### Based on the initial predicted labels, reconsturct the reference and update labels:
```
python -m pipelines.TORC_pipeline --torc_step TORC_reconstruct --select_method $fs_method --res_dir $RES_DIR --data_dir $DATA_DIR --ref_dir $REF_DIR --NMLP_dir $NMLP_DIR --sample_method $smethod --expand_ref $do_expand --Entropy_dir $EP_DIR --scale
```

### Parameters

The TORC_pipeline script provides several configurable command-line arguments to customize its execution. Below is a description of each argument and its possible values.

- '--torc_step' (str, required)
  - Specifies the step of the TORC pipeline.
  - Choices: "TORC_initial", "TORC_reconstruct".
- '--select_method' (str, required)
  - Defines the feature selection method.
  - Choices: "noFS" (no feature selection), "F-test", "Seurat".
- '--n_features' (int, optional, default: 1000)
  - Number of features selected..
- '--res_dir'(str, required)
  - Directory where results will be stored.
- '--ref_dir' (str, required)
  - Directory containing the reference dataset.
- '--data_dir' (str, required)
  - Directory containing the target dataset.
- '-train_sample' (str, optional)
  - Train sample ID for downsampling
- '-test_sample' (str, optional)
  - Target sample ID for downsampling
- '--NMLP_dir' (str, optional)
  - Directory containing Initial MLP results.
- '--Entropy_dir' (str, optional)
  - Directory containing entropy status results.
- '--expand_ref' (bool, optional)
  - Specifies whether to expand the reference dataset.
- '--sample_method', help="Define Sampling method")
  - Defines the sampling method for reference reconsturction to be used.
  - Choices: "Prop" (proportion-based sampling), "Prob" (probability-based sampling)
- '--sample_size' (int, optional, default: 1000)
  - Approximate sample size for analysis.
- '--rseed' (int, optional, default: 1234)
  - Specifies the random seed for reproducibility.
- '--scale' (flag argument)
  - Enables data scaling if specified.   







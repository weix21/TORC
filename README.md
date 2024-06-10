# TORC:Target-Oriented Reference Construction for supervised cell-type identification in scRNA-seq

TORC aims to construct a reference suitable for a given target data. It first uses an off-the-shelf supervised method to label the target cells, from which TORC first estimates cell-type composition in the target. TORC then add target cells with high-confidence labels to the reference to form an expanded reference pool. Then TORC resamples from the pool to construct a new reference according to the target composition. The reconstructed reference is used to build the final classifier. 

![Figure1](https://github.com/weix21/TORC/assets/81911479/38dcfd85-7f27-4159-9778-2f9fdf21718a)

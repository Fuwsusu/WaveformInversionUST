# FuWs Utilities

This folder contains custom workflow scripts and helpers, including sparse-array tools for Ali in-vivo breast datasets.

## Sparse-array helpers in this folder

- `downsampleAliInvivoElements.m`  
  File-based workflow. Loads a `.mat` dataset and writes a sparse `.mat` dataset.
- `makeSparseRingArrayData.m`  
  In-memory callable function. Returns sparse `full_dataset`, sparse `transducerPositionsXY`, and `keepIdx`.

## Quick start

```matlab
[out, keepIdx] = downsampleAliInvivoElements('SampleData/Malignancy.mat', ...
    'SampleData/Malignancy_sparse256.mat', 256);
```

```matlab
[full_dataset_sparse, transducerPositionsXY_sparse, keepIdx] = ...
    makeSparseRingArrayData(full_dataset, transducerPositionsXY, 256);
```

See `FuWs/SampleData/README.md` for folder-specific notes.

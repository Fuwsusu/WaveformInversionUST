# FuWs/SampleData

This subfolder is created for FuWs-local documentation and sample-data notes.

If you use the sparse-array utilities in `FuWs/`:

- Input files are typically from repository root `SampleData/`.
- Output sparse files can be saved either to repository root `SampleData/` or inside this folder.

Example:

```matlab
[out, keepIdx] = downsampleAliInvivoElements('SampleData/BenignCyst.mat', ...
    'FuWs/SampleData/BenignCyst_sparse256.mat', 256);
```

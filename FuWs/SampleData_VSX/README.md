# FuWs/SampleData_VSX

Utilities for building sparse-array datasets from VSX phantom MAT files.

## Files

- `downsampleVSXElements_fws.m`: callable function for VSX sparse-element conversion.
- `run_downsampleVSXElements_fws.m`: run script with the VSX phantom path provided by user.

## Quick start

In MATLAB:

```matlab
run('FuWs/SampleData_VSX/run_downsampleVSXElements_fws.m')
```

This will read:

- `D:\Document_ING_fws\WaveformInversionUST\SampleData\VSX_YezitronixPhantom1.mat`

and write output into:

- `FuWs/SampleData_VSX/VSX_YezitronixPhantom1_sparse256.mat`

function [outputStruct, keepIdx] = downsampleVSXElements_fws(inputFile, outputFile, targetNumElements)
% downsampleVSXElements_fws
% Build sparse-ring channel data from a VSX phantom MAT file.
%
% Usage:
%   [out, keepIdx] = downsampleVSXElements_fws(...
%       'SampleData/VSX_YezitronixPhantom1.mat', ...
%       'FuWs/SampleData_VSX/VSX_YezitronixPhantom1_sparse256.mat', 256);
%
% Inputs:
%   inputFile          - Source VSX MAT-file path
%   outputFile         - Output sparse MAT-file path
%   targetNumElements  - Number of ring elements to keep
%
% Outputs:
%   outputStruct       - Struct saved to outputFile
%   keepIdx            - Kept element indices
%
% Notes:
%   This utility expects VSX data fields used by MultiFrequencyWaveformInvVSX.m:
%   - full_dataset: [Nt x Nrx x Ntx]
%   - transducerPositionsXY: [2 x N]
%   - time: [Nt x 1] or [1 x Nt] (recommended)

    arguments
        inputFile (1,:) char
        outputFile (1,:) char
        targetNumElements (1,1) double {mustBeInteger,mustBePositive}
    end

    in = load(inputFile);

    requiredFields = {'full_dataset', 'transducerPositionsXY'};
    for k = 1:numel(requiredFields)
        if ~isfield(in, requiredFields{k})
            error('Missing required field "%s" in %s.', requiredFields{k}, inputFile);
        end
    end

    if ~isfield(in, 'time')
        warning('Field "time" is missing in %s. The output will still be written.', inputFile);
    end

    originalNumElements = size(in.transducerPositionsXY, 2);

    helperDir = fullfile(fileparts(mfilename('fullpath')), '..', 'SampleData');
    if exist('makeSparseRingArrayData_fws', 'file') ~= 2
        addpath(helperDir);
    end

    [in.full_dataset, in.transducerPositionsXY, keepIdx] = makeSparseRingArrayData_fws(...
        in.full_dataset, in.transducerPositionsXY, targetNumElements);

    in.sparseArray.originalNumElements = originalNumElements;
    in.sparseArray.targetNumElements = targetNumElements;
    in.sparseArray.keepIdx = keepIdx;
    in.sparseArray.createdBy = mfilename;
    in.sparseArray.createdAtUTC = char(datetime('now', 'TimeZone', 'UTC', ...
        'Format', 'yyyy-MM-dd HH:mm:ss z'));

    outDir = fileparts(outputFile);
    if ~isempty(outDir) && ~exist(outDir, 'dir')
        mkdir(outDir);
    end

    save(outputFile, '-struct', 'in', '-v7.3');

    fprintf('Saved sparse VSX dataset: %s\n', outputFile);
    fprintf('Element count: %d -> %d\n', originalNumElements, targetNumElements);

    outputStruct = in;
end

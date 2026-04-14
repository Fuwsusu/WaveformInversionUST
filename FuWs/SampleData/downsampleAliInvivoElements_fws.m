function [outputStruct, keepIdx] = downsampleAliInvivoElements_fws(inputFile, outputFile, targetNumElements)
% downsampleAliInvivoElements
% Reduce Ali's in-vivo KCI breast dataset to a sparse ring-array subset.
%
% Usage:
%   [out, keepIdx] = downsampleAliInvivoElements('SampleData/Malignancy.mat', ...
%       'SampleData/Malignancy_sparse256.mat', 256)
%
% Inputs:
%   inputFile          - Path to source MAT file (e.g., BenignCyst/Malignancy)
%   outputFile         - Path for sparse output MAT file
%   targetNumElements  - Number of ring elements to keep (e.g., 512/256/128)
%
% Outputs:
%   outputStruct       - MAT struct that was written to outputFile
%   keepIdx            - Kept ring element indices
%
% Notes:
%   1) This function assumes data fields are compatible with
%      MultiFrequencyWaveformInvKCI.m:
%      - full_dataset: [Nt x Nrx x Ntx]
%      - transducerPositionsXY: [2 x N]
%      - time: [Nt x 1] or [1 x Nt]
%   2) Uses callable helper: FuWs/makeSparseRingArrayData.m

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

    originalNumElements = size(in.transducerPositionsXY, 2);


    if exist('makeSparseRingArrayData', 'file') ~= 2
        thisDir = fileparts(mfilename('fullpath'));
        addpath(thisDir);
    end

    [in.full_dataset, in.transducerPositionsXY, keepIdx] = makeSparseRingArrayData(...
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

    fprintf('Saved sparse dataset: %s\n', outputFile);
    fprintf('Element count: %d -> %d\n', originalNumElements, targetNumElements);

    outputStruct = in;
end

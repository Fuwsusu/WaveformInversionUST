function [full_dataset_sparse, transducerPositionsXY_sparse, keepIdx] = ...
    makeSparseRingArrayData(full_dataset, transducerPositionsXY, targetNumElements)
% makeSparseRingArrayData
% Callable helper to sparsify ring-array UST channel data in memory.
%
% Usage:
%   [full_dataset_sparse, transducerPositionsXY_sparse, keepIdx] = ...
%       makeSparseRingArrayData(full_dataset, transducerPositionsXY, 256);
%
% Inputs:
%   full_dataset           - [Nt x Nrx x Ntx] full ring-array channel data
%   transducerPositionsXY  - [2 x N] transducer coordinates in meters
%   targetNumElements      - Number of elements to retain on the ring
%
% Outputs:
%   full_dataset_sparse        - [Nt x N_sparse x N_sparse]
%   transducerPositionsXY_sparse - [2 x N_sparse]
%   keepIdx                    - Indices of retained ring elements

    arguments
        full_dataset (:,:,:) {mustBeNumeric}
        transducerPositionsXY (2,:) {mustBeNumeric}
        targetNumElements (1,1) double {mustBeInteger,mustBePositive}
    end

    numElements = size(transducerPositionsXY, 2);

    if size(full_dataset, 2) ~= numElements || size(full_dataset, 3) ~= numElements
        error(['Expected full_dataset size [Nt x N x N] with N matching transducerPositionsXY. ', ...
               'Got size(full_dataset) = [%s], N = %d.'], num2str(size(full_dataset)), numElements);
    end

    if targetNumElements > numElements
        error('targetNumElements (%d) cannot exceed original element count (%d).', ...
            targetNumElements, numElements);
    end

    keepIdx = unique(round(linspace(1, numElements, targetNumElements)));
    if numel(keepIdx) < targetNumElements
        step = floor(numElements / targetNumElements);
        keepIdx = 1:step:numElements;
        keepIdx = keepIdx(1:targetNumElements);
    end

    full_dataset_sparse = full_dataset(:, keepIdx, keepIdx);
    transducerPositionsXY_sparse = transducerPositionsXY(:, keepIdx);
end

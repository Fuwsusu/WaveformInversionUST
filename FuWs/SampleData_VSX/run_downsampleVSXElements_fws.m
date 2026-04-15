% run_downsampleVSXElements_fws
% 针对 VSX 真实体模数据构建稀疏阵元数据。
%
% 你提供的输入路径：
%   D:\Document_ING_fws\WaveformInversionUST\SampleData\VSX_YezitronixPhantom1.mat
%
% 建议：
% 1) 在 MATLAB 中直接运行本脚本。
% 2) 按需修改 targetNumElements（例如 512 / 256 / 128）。

clear; clc;

% ====== 用户参数 ======
inputFile = 'D:\Document_ING_fws\WaveformInversionUST\SampleData\VSX_YezitronixPhantom1.mat';
targetNumElements = 256;

% 在仓库中新建输出目录 FuWs/SampleData_VSX
repoRoot = fileparts(fileparts(fileparts(mfilename('fullpath'))));
outputDir = fullfile(repoRoot, 'FuWs', 'SampleData_VSX');
if ~exist(outputDir, 'dir')
    mkdir(outputDir);
end

[~, inputName, ~] = fileparts(inputFile);
outputFile = fullfile(outputDir, sprintf('%s_sparse%d.mat', inputName, targetNumElements));

% ====== 调用函数 ======
[outputStruct, keepIdx] = downsampleVSXElements_fws(inputFile, outputFile, targetNumElements);

% ====== 简要结果输出 ======
fprintf('\nDone.\n');
fprintf('Input : %s\n', inputFile);
fprintf('Output: %s\n', outputFile);
fprintf('Kept elements: %d\n', numel(keepIdx));

if isfield(outputStruct, 'sparseArray')
    disp(outputStruct.sparseArray);
end

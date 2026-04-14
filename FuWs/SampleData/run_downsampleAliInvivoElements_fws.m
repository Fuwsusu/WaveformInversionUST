% run_downsampleAliInvivoElements_fws
% 调用 downsampleAliInvivoElements_fws 对乳腺数据做稀疏阵列降采样。
%
% 使用说明：
% 1) 直接运行本脚本。
% 2) 根据需要修改 targetNumElements（如 512/256/128）。
%
% 输入数据（你提供的路径）：
%   D:\Document_ING_fws\WaveformInversionUST\SampleData\BenignCyst.mat

clear; clc;

% ====== 用户参数 ======
inputFile = 'D:\Document_ING_fws\WaveformInversionUST\SampleData\BenignCyst.mat';
targetNumElements = 256;

% 输出文件：自动在同目录下生成 *_sparseXXX.mat
[inputDir, inputName, ~] = fileparts(inputFile);
outputFile = fullfile(inputDir, sprintf('%s_sparse%d.mat', inputName, targetNumElements));

% ====== 调用函数 ======
[outputStruct, keepIdx] = downsampleAliInvivoElements_fws(inputFile, outputFile, targetNumElements);

% ====== 简要结果输出 ======
fprintf('\nDone.\n');
fprintf('Input : %s\n', inputFile);
fprintf('Output: %s\n', outputFile);
fprintf('Kept elements: %d\n', numel(keepIdx));

% 可选：检查输出字段
if isfield(outputStruct, 'sparseArray')
    disp(outputStruct.sparseArray);
end

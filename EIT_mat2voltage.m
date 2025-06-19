function EIT_mat2voltage()
% EIT_MAT2VOLTAGE  Batch convert strain MAT → voltage & diff‑recon  (no inputs)
% -------------------------------------------------------------------------
% 固定参数：如需调整直接改下面常量
% -------------------------------------------------------------------------
DIR_PATH     = "/home/kairui/projects/EIT/SlideTraj";                % *.mat 目录
MODEL_TAG    = 'h2c';        % EIDORS 模型标签（16‑electrode circular）
GF           = 8;            % Gauge factor
SIGMA0       = 1;            % 背景电导率
BASELINE_IDX = 1;            % 差分重建基准帧
SAVE_PNG     = false;         % 是否保存 *_recon.png
USE_PARFOR   = false;        % 是否并行文件循环
% -------------------------------------------------------------------------

% ---------- 1. 启动 EIDORS ----------
eidors_start_file = '~/projects/EIDORs/eidors-v3.12-ng/eidors/startup.m';
if exist(eidors_start_file,'file')
    run(eidors_start_file);
else
    error("EIDORS startup file not found:\n  %s",eidors_start_file);
end
inv_model = buildEidorsDiffModel(16,MODEL_TAG);      % 相邻激励差分模型
fprintf("EIDORS initialised.\n");

% ---------- 2. 扫描 MAT 文件 ----------
files = dir(fullfile(DIR_PATH,"*.mat"));
if isempty(files)
    error("No MAT files found in %s",DIR_PATH);
end
fprintf("Found %d MAT files in %s\n",numel(files),DIR_PATH);

% ---------- 3. 循环处理 ----------
loopFn = @(idx) processOne(files(idx),MODEL_TAG,GF,SIGMA0,...
                           BASELINE_IDX,inv_model,SAVE_PNG);

if USE_PARFOR
    parfor i = 1:numel(files), loopFn(i); end
else
    for i = 1:numel(files),  loopFn(i); end
end
fprintf("\n✓ All files processed.\n");
end
%% ========================================================================
function processOne(finfo,model_tag,GF,sigma0,...
                    baselineIdx,inv_model,savePNG)
srcMat = fullfile(finfo.folder,finfo.name);
S      = load(srcMat,"strain","tlist");

if ~isfield(S,"strain")
    warning("Skip %s (no strain variable)",finfo.name); return;
end
bigArr  = S.strain;                 % [T,N,3]
timeVec = parseTimeVector(S);

fprintf("→ %-20s  (%d frames)\n",finfo.name,size(bigArr,1));

% ---- 计算电压 -----------------------------------------------------------
[~, allV] = computeSigmaAndVoltage(bigArr,GF,sigma0,model_tag);

dstVol = fullfile(finfo.folder,replace(finfo.name,".mat","_voltage.mat"));
save(dstVol,"allV","timeVec");              % 保存电压
fprintf("   voltage saved -> %s\n",dstVol);

% ---- 差分重建 -----------------------------------------------------------
doDiffReconstructionFromVoltage(allV,baselineIdx,inv_model,timeVec);

if savePNG
    pngFile = fullfile(finfo.folder,replace(finfo.name,".mat","_recon.png"));
    saveas(gcf,pngFile);                    % 保存图像
    fprintf("   recon saved   -> %s\n",pngFile);
end
close(gcf);                                 % 关闭窗口
end
%% ========================================================================
% 其余函数（computeSigmaAndVoltage、buildEidorsDiffModel、
% doDiffReconstructionFromVoltage、parseTimeVector）保持不变，
% 请确保它们位于同一文件夹或已加入 MATLAB path。


%% ------------------------------------------------------------------------
function [bigArr, N, t] = readStrainCsvData(csv_file)
% [bigArr, N, t] = readStrainCsvData(csv_file)
%   读取CSV, 返回 bigArr(t,N,3)
%   对应: x,y坐标, 以及 第i步的应变

    fid = fopen(csv_file,'r');
    if fid<0
        error('could not open the file: %s', csv_file);
    end

    % 跳过前4行
    for skip=1:4
        fgetl(fid);
    end

    % 第5行 => "Nodes, N"
    line = fgetl(fid);
    tokens = regexp(line, '^%?\s*Nodes,\s*(\d+)', 'tokens');
    if isempty(tokens)
        error('could not load the 5th node: %s', line);
    end
    N = str2double(tokens{1}{1});

    % 第6行 => "Expressions, t"
    line = fgetl(fid);
    tokens = regexp(line, '^%?\s*Expressions,\s*(\d+)', 'tokens');
    if isempty(tokens)
        error('无法在第6行解析时间步: %s', line);
    end
    t = str2double(tokens{1}{1});

    % 跳过第7,8行
    fgetl(fid); fgetl(fid);
    % 跳过第9行(列标题)
    fgetl(fid);

    % 读取N行数据 => [x,y,z, val_1..val_t]
    data = zeros(N, 3+t);
    for n=1:N
        line = fgetl(fid);
        if ~ischar(line)
            error('文件行数不足, 期望N=%d行, 实际只到%d行.', N, n-1);
        end
        vals = str2num(line); %#ok<ST2NM>
        if length(vals)<3+t
            error('第%d行列数不足(需要 %d, 实际 %d).', n, 3+t, length(vals));
        end
        data(n,:) = vals(1:(3+t));
    end
    fclose(fid);

    % 构造 bigArr => [t, N, 3]
    x = data(:,1);
    y = data(:,2);
    valMatrix = data(:,4:end);  % [N x t]

    bigArr = zeros(t, N, 3);
    for i=1:t
        bigArr(i,:,1) = x';               % x
        bigArr(i,:,2) = y';               % y
        bigArr(i,:,3) = valMatrix(:,i)';  % 第 i 列(第 i 个时间步)
    end

    fprintf('读取完毕: N=%d, t=%d.\n', N, t);
    fprintf('bigArr大小= [%d x %d x 3].\n', size(bigArr,1), size(bigArr,2));
end


%% ------------------------------------------------------------------------
function [sigma_elem_all, allVoltage, globalMin, globalMax, fmdl] = ...
    computeSigmaAndVoltage(bigArr, GF, sigma0, model_tag)
% 计算电导率(插值)并做前向求解, 返回:
%   sigma_elem_all: cell(t,1), 每步单元电导率
%   allVoltage    : cell(t,1), 每步电极电压测量
%   globalMin, globalMax: 全局电导率最小/最大
%   fmdl          : EIDORS前向模型

    [t, N, ~] = size(bigArr);

    % 1) 创建EIDORS模型
    imdl = mk_common_model(model_tag,16);
    fmdl = imdl.fwd_model;
    fmdl.output = 'all';  % 若版本支持输出节点电位

    nodes = fmdl.nodes;
    elems = fmdl.elems;

    allVoltage = cell(t,1);
    sigma_elem_all = cell(t,1);

    globalMin = inf;
    globalMax = -inf;

    for i=1:t
        x_i = bigArr(i,:,1)';  % [N x 1]
        y_i = bigArr(i,:,2)';
        strain_i = bigArr(i,:,3)';

        % 应变 => 电导率
        sigma_i = sigma0.*(1 - GF.*strain_i);

        % 散点插值 => 单元
        F = scatteredInterpolant(x_i, y_i, sigma_i, 'linear','nearest');
        elem_centers = (nodes(elems(:,1),:)+nodes(elems(:,2),:)+nodes(elems(:,3),:))/3;
        sigma_elem = F(elem_centers(:,1), elem_centers(:,2));

        % 更新全局 min/max
        localMin = min(sigma_elem);
        localMax = max(sigma_elem);
        if localMin<globalMin, globalMin=localMin; end
        if localMax>globalMax, globalMax=localMax; end

        sigma_elem_all{i} = sigma_elem;

        % 前向求解 => 电压
        img = mk_image(fmdl, sigma_elem);
        vh = fwd_solve(img);
        allVoltage{i} = vh.meas;

        fprintf('Step %d/%d: %d通道\n', i, t, length(vh.meas));
    end
end


%% ------------------------------------------------------------------------
function plotSigmaDistribution(sigma_elem_all, fmdl, globalMin, globalMax)
% plotSigmaDistribution
%   在同一个figure里, subplot方式绘制每个时间步电导率分布

    t = length(sigma_elem_all);
    figure('Name','Sigma distribution','NumberTitle','off');
    for i=1:t
        subplot(ceil(sqrt(t)), ceil(sqrt(t)), i);

        img = mk_image(fmdl, sigma_elem_all{i});
        show_fem(img,[1,0, globalMin, globalMax]);
        % colorbar;  %若需要
        title(['Conductivity, step=' num2str(i)]);
    end
end


%% ------------------------------------------------------------------------
function plotVoltageGroups(allVoltage, timeVec)
% plotVoltageGroups
%   将电压(通道x时间)做分箱, 并在同一figure subplot里绘制.
%   1) 构造 volMatrix(t x M)
%   2) 求channelAvg
%   3) 分箱 => subplot

    tCount = length(allVoltage);
    M = length(allVoltage{1});

    % 构造 volMatrix
    volMatrix = zeros(tCount,M);
    for i=1:tCount
        volMatrix(i,:) = allVoltage{i}';
    end

    % 计算平均值
    channelAvg = mean(volMatrix,1);
    minVal = min(channelAvg);
    maxVal = max(channelAvg);

    binWidth=0.05;  % 例如 0.05 V
    nBins = ceil((maxVal-minVal)/binWidth);

    groups = cell(nBins,1);
    for c=1:M
        offset = channelAvg(c)-minVal;
        binIndex = floor(offset/binWidth)+1;
        binIndex = max(binIndex,1);
        binIndex = min(binIndex,nBins);
        groups{binIndex} = [groups{binIndex}, c];
    end

    figure('Name','Voltage grouping','NumberTitle','off');
    rows = ceil(sqrt(nBins));
    cols = ceil(nBins/rows);

    for b=1:nBins
        subplot(rows,cols,b);
        chInBin = groups{b};
        if isempty(chInBin)
            text(0.5,0.5,'No channels','HorizontalAlignment','center');
            axis off;
        else
            plot(timeVec, volMatrix(:,chInBin), 'LineWidth',1.2);
            grid on;
            lowB = minVal + (b-1)*binWidth;
            highB= minVal + b*binWidth;
            title(sprintf('Bin %d: [%.2f, %.2f)', b, lowB, highB));
        end
    end
    fprintf('plotVoltageGroups: %d bins, each subplot one bin.\n', nBins);
end
%% ------------------------------------------------------------------------
function inv_model = buildEidorsDiffModel(n_elec,model_tag)
% buildEidorsDiffModel: 构建相邻激励+差分 EIDORS 模型
    fmdl_struct = mk_common_model(model_tag, n_elec);
    fmdl = fmdl_struct.fwd_model;
    fmdl.type = 'fwd_model';

    % 相邻激励+差分
    stim = mk_stim_patterns(n_elec, 1, [0,1],[0,1],{'do_redundant'},1);
    % scale if needed
    fmdl.stimulation = stim;

    inv_model = eidors_obj('inv_model','EIT_diff_imaging');
    inv_model.fwd_model = fmdl;
    inv_model.reconst_type = 'difference';
    inv_model.solve = @inv_solve_diff_pdipm; 
    inv_model.hyperparameter.value = 1e-2;
    inv_model.RtR_prior = @prior_laplace;
    inv_model.meas_constraint = @meas_icov;

    % 背景图像
    img_bkg = mk_image(fmdl,1);
    inv_model.jacobian_bkgnd = img_bkg;
end
%% ------------------------------------------------------------------------
function img_rec_cell = doDiffReconstructionFromVoltage(allVoltage, baselineIdx, inv_model, timeVec)
% doDiffReconstructionFromVoltage
%   使用 EIDORS 差分成像( inv_model )对 allVoltage 做差分重建.
%   在同一figure中以 subplot 形式显示每个时刻的重建结果.
%
% 输入:
%   allVoltage  : cell(t,1), 每个时间步的电极测量向量 (M x 1)
%   baselineIdx : 基准帧索引(1..t)
%   inv_model   : EIDORS 逆问题模型 (reconst_type='difference')
%   timeVec     : (t x 1)时间向量, 用于给图加标题(可为空)
%
% 输出:
%   img_rec_cell: cell(t,1), 每个时刻的重建图像对象(EIDORS image)
%
% 使用:
%   load('allVoltage.mat','allVoltage');
%   baselineIdx=1;
%   timeVec=linspace(0,1.5,length(allVoltage))';
%   inv_model=... (差分模型)
%   doDiffReconstructionFromVoltage(allVoltage, baselineIdx, inv_model, timeVec);

    tCount = length(allVoltage);
    if baselineIdx<1 || baselineIdx>tCount
        error('baselineIdx=%d 超出时间步范围1..%d', baselineIdx, tCount);
    end

    % 1) 取baseline帧电压
    v0 = allVoltage{baselineIdx};  % (M x 1)

    % 2) 预分配存储重建图像
    img_rec_cell = cell(tCount,1);

    % 3) 创建 figure
    figure('Name','Diff Reconstruction from allVoltage','NumberTitle','off');
    for i = 1:tCount
        % 3.1) 取第 i 帧电压
        v_t = allVoltage{i};

        % 3.2) 做差分重建
        %  EIDORS 差分 => inv_solve(inv_model, v0, v_t)
        img_rec = inv_solve(inv_model, v0, v_t);
        img_rec_cell{i} = img_rec;

        % 3.3) subplot绘图
        subplot(ceil(sqrt(tCount)), ceil(sqrt(tCount)), i);
        show_fem(img_rec, [1,0]);  % 只显示颜色(不显示网格)
        if nargin>=4 && ~isempty(timeVec) && length(timeVec)==tCount
            title(sprintf('t=%.3g s', timeVec(i)));
        else
            title(sprintf('frame %d', i));
        end
    end

    fprintf('已完成差分重建, 并在同一figure用 subplot 显示 %d 帧.\n', tCount);
end

function tvec = parseTimeVector(S)
% If tlist exists (e.g. 'range(0,0.03,3)'), build a numeric vector.
if isfield(S,'tlist') && ischar(S.tlist) && startsWith(S.tlist,'range')
    num = sscanf(S.tlist,'range(%f,%f,%f)');
    tvec = num(1):num(2):num(3);
    tvec = tvec(:);
else
    % fallback: 0,1,2,…
    T = size(S.strain,1);
    tvec = (0:T-1).';
end
end

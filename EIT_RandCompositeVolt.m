% ========================================================================
%  EIT_RandCompositeVolt.m   (updated 2025‑05‑07)
%  ---------------------------------------------------------
%  • 随机组合  Slide  +  单点动作(Press / Tap / Punch)    →    电压合并
%  • Slide 现支持三类轨迹: 直线(line) / 正弦(sine) / 圆(circle)
%  • 每个动作单独求解 (各自网格)  →  导出 Cut‑Plane Strain CSV
%  • Strain  →  Conductivity  →  Electrode Voltage (EIDORS 前向)
%  • 把各动作的电压序列、时间轴拼接保存到一个 MAT
% ========================================================================
%  依赖：
%  - LiveLink for MATLAB (COMSOL) 已配置  (javaaddpath 已在安装脚本处理)
%  - EIDORS 已安装，修改下方 eidors_startup_path
% ------------------------------------------------------------------------

%% =====  环境路径  =====
addpath('/mnt/disk2/software/comsol/mli');           % LiveLink   (★修改)
eidors_startup_path = '/home/kairui/projects/EIDORs/eidors-v3.12-ng/eidors/startup.m';
run(eidors_startup_path);                            % 启动 EIDORS
import com.comsol.model.util.ModelUtil               % COMSOL LiveLink API

%% =====  用户参数  =====
base_single = '/home/kairui/projects/EIT/Comsol6.2/EIT3D_SingleMotion.mph';
base_slide  = '/home/kairui/projects/EIT/Comsol6.2/EIT3D_Slide.mph';
out_dir     = '/home/kairui/projects/EIT/RandCompositeOut';
if ~exist(out_dir,'dir'), mkdir(out_dir); end

SEQ_LEN         = 4;                                % 每条序列动作数
MAX_PER_KIND    = 2;                                % 同类型动作出现上限
SINGLE_TYPES    = {'Press','Tap','Punch'};          % 单点动作集合

% ---- Slide 轨迹参数预设 -------------------------------------------------
%  直线   : x(t)=x0+k_x*t,  y(t)=y0+k_y*t
%  正弦   : 轴向周期漂移  (axis = 'x' 或 'y')
%  圆形   : 半径 R, 角速度 w (rad/s)
SLIDE_PARAM_SET = {
    % ----- 直线(Line)  -----
    struct('type','line',  'k_x', 0.001,'k_y',  0.000);
    struct('type','line',  'k_x', 0.001,'k_y', -0.001); 
    struct('type','line',  'k_x', 0.000,'k_y',  0.001); 
    % ----- 正弦(Sine)  -----
    struct('type','sine',  'axis','x', 'A', 0.003,'w', 2*pi/3);
    struct('type','sine',  'axis','y', 'A', 0.003,'w', 2*pi/3);
    % ----- 圆形(Circle) -----
    struct('type','circle','R',  0.003,'w', 2*pi/3)
};

% 位置与力学常量
x0_range   = [-0.005, 0.005];   % m
y0_range   = [-0.005, 0.005];   % m
a_val      = 100;               % exp 衰减
force_val  = 1e9;               % 载荷幅值

tlist_single = 'range(0,0.1,3)';
tlist_slide  = 'range(0,0.1,3)';

% Strain → Voltage 参数
GF        = 8;
SIGMA0    = 1;
MODEL_TAG = 'h2C';              % 4096‑elem, 16‑elec CEM 2‑D 圆模型
rng('shuffle');                 % 随机种子

%% =====  随机生成动作序列  =====
x0 = randRange(x0_range);
y0 = randRange(y0_range);

motions = {};  cntSlide = 0;     cntSingle = containers.Map(SINGLE_TYPES,zeros(1,numel(SINGLE_TYPES)));

while numel(motions) < SEQ_LEN
    if rand < .4 && cntSlide < MAX_PER_KIND   % 尝试放置 Slide
        cntSlide = cntSlide + 1;
        mv.kind  = 'Slide';
        mv.param = SLIDE_PARAM_SET{ randi(numel(SLIDE_PARAM_SET)) };
        motions{end+1} = mv; %#ok<SAGROW>
    else                                       % 尝试放置单点动作
        st = SINGLE_TYPES{ randi(numel(SINGLE_TYPES)) };
        if cntSingle(st) < MAX_PER_KIND
            cntSingle(st) = cntSingle(st)+1;
            mv.kind  = st;
            mv.param = [];
            motions{end+1} = mv; %#ok<SAGROW>
        end
    end
end

fprintf('◆ Random Sequence:\n');
for i = 1:numel(motions)
    m = motions{i};
    if strcmp(m.kind,'Slide')
        fprintf('   - Slide (%s)\n', m.param.type);
    else
        fprintf('   - %s\n', m.kind);
    end
end

%% =====  预生成 EIDORS 前向模型 (共享) =====
tmp = mk_common_model(MODEL_TAG,16);
fmdl_global = tmp.fwd_model;

%% =====  全局容器初始化 =====
all_voltage = {};      % 每帧电压  cell{T×1}
all_time    = [];      % 全局时间向量
motion_tag  = {};      % 每帧标签

%% =====  循环求解各动作 =====
for idx = 1:numel(motions)
    mov = motions{idx};
    mk  = mov.kind;      % 'Slide' / 'Press' / ...

    % -------- 载入相应模板 mph --------
    if strcmp(mk,'Slide')
        model = mphload(base_slide);
    else
        model = mphload(base_single);
    end

    % -------- 统一参数 --------
    model.param.set('x0', num2str(x0));
    model.param.set('y0', num2str(y0));
    model.param.set('a',  num2str(a_val));
    model.param.set('force_scale', num2str(force_val));

    % -------- 各动作专属设置 --------
    if strcmp(mk,'Slide')
        sp = mov.param;
        switch sp.type
            case 'line'    % 直线: x = x0 + k_x*t, y = y0 + k_y*t
                model.func('an1').set('expr', sprintf('%g+%g*t', x0, sp.k_x));
                model.func('an2').set('expr', sprintf('%g+%g*t', y0, sp.k_y));
            case 'sine'    % 正弦: 单轴振荡
                if strcmp(sp.axis,'x')
                    model.func('an1').set('expr', sprintf('%g+%g*sin(%g*t)', x0, sp.A, sp.w));
                    model.func('an2').set('expr', sprintf('%g', y0));
                else % 'y'
                    model.func('an1').set('expr', sprintf('%g', x0));
                    model.func('an2').set('expr', sprintf('%g+%g*sin(%g*t)', y0, sp.A, sp.w));
                end
            case 'circle'  % 圆形: 半径 R
                model.func('an1').set('expr', sprintf('%g+%g*cos(%g*t)', x0, sp.R, sp.w));
                model.func('an2').set('expr', sprintf('%g+%g*sin(%g*t)', y0, sp.R, sp.w));
            otherwise
                error('Unknown slide trajectory type: %s', sp.type);
        end
        model.study('std1').feature('time').set('tlist', tlist_slide);
        tvec = evalTlist(tlist_slide);
    else % -------- 单点动作 --------
        ftag = ['F', lower(mk)];   % 'Fpress' / 'Ftap' / 'Fpunch'
        expr = sprintf('- %s(t)*exp(-a*((x-x0)^2+(y-y0)^2))*%.0e', ftag, force_val);
        model.physics('solid').feature('bndl1') ...
             .set('FperArea',{expr,'0','0'});
        model.study('std1').feature('time').set('tlist', tlist_single);
        tvec = evalTlist(tlist_single);
    end
    model.study('std1').feature('time').set('useinitsol','off');

    % -------- 求解 COMSOL --------
    model.study('std1').run;
    fprintf('     ✔ %s solved (%d frames)\n', mk, numel(tvec));

    % -------- 导出 Strain CSV --------
    tmp_csv = fullfile(out_dir,'tmp.csv');
    if any(string(model.result.export.tags)=="expCsv")
        model.result.export.remove('expCsv');
    end
    ex = model.result.export.create('expCsv','Data');
    ex.set('data','cpl1');
    ex.set('expr',{'solid.ep1'});
    ex.set('filename',tmp_csv);  ex.run;
    fprintf('     ✔ exported\n');
    strain = readCsv(tmp_csv); delete(tmp_csv);

    % -------- Strain → Voltage (EIDORS) --------
    [~, Vcell] = computeSigmaAndVoltage(strain, GF, SIGMA0, fmdl_global);

    % -------- 合并电压和时间线 --------
    if isempty(all_time)
        offset = 0;
    else
        offset = all_time(end);
    end
    all_time    = [all_time;  tvec(:)+offset];          %#ok<AGROW>
    all_voltage = [all_voltage; Vcell];                 %#ok<AGROW>
    motion_tag  = [motion_tag; repmat({mk}, numel(Vcell),1)]; %#ok<AGROW>

    % -------- 释放模型 --------
    ModelUtil.remove('model');
    java.lang.System.gc(); pause(0.1);
end

%% =====  保存结果 =====
case_id   = datestr(now,'yyyymmdd_HHMMSS');
save_file = fullfile(out_dir, sprintf('Composite_%s.mat',case_id));
save(save_file,'all_voltage','all_time','motion_tag',...
               'x0','y0','a_val','force_val','motions');
fprintf('\n★ Saved voltage MAT → %s\n', save_file);

%% =====  可选：差分重建可视化 =====
INV_MODEL = buildEidorsDiffModel(16);
doDiffReconstructionFromVoltage(all_voltage,1,INV_MODEL,all_time);

%% ==========  分箱绘制电压曲线(同一个figure, subplot) ==========
%  你可以定义时间向量, e.g. 0~12 s
    timeVec = linspace(0,12, length(all_voltage))';
    plotVoltageGroups(all_voltage, timeVec);
% ========================================================================
%                            HELPER FUNCTIONS
% ========================================================================

function v = randRange(r)
% uniform random in [r(1), r(2)]
v = r(1) + (r(2)-r(1))*rand;
end
% ------------------------------------------------------------------------
function t = evalTlist(str)
% convert 'range(a,step,b)' to numeric column vector
num = sscanf(str,'range(%f,%f,%f)');
t   = (num(1):num(2):num(3)).';
end
% ------------------------------------------------------------------------
function strain = readCsv(csvFile)
% Read COMSOL Cut‑Plane CSV → strain(T,N,3)
fid=fopen(csvFile,'r'); for i=1:4, fgetl(fid); end
N = sscanf(fgetl(fid),'%% Nodes,%d');
T = sscanf(fgetl(fid),'%% Expressions,%d');
for i=1:3, fgetl(fid); end
raw=zeros(N,3+T);
for n=1:N, raw(n,:)=sscanf(fgetl(fid),'%f,').'; end
fclose(fid);
x=raw(:,1); y=raw(:,2); vs=raw(:,4:end);
strain=zeros(T,N,3);
for t=1:T
    strain(t,:,1)=x.'; strain(t,:,2)=y.'; strain(t,:,3)=vs(:,t).';
end
end

% ------------------------------------------------------------------------
function [sigma_elem_all, allV] = computeSigmaAndVoltage(bigArr, GF, sigma0, fmdl)
% bigArr: [T,N,3]; return allV{t} each (M×1)
[T,~,~]=size(bigArr);
nodes=fmdl.nodes; elems=fmdl.elems;
sigma_elem_all=cell(T,1);   allV=cell(T,1);
for i=1:T
    xi=bigArr(i,:,1)'; yi=bigArr(i,:,2)'; eps=bigArr(i,:,3)';
    sigma_n = sigma0 .* (1 - GF .* eps);
    F = scatteredInterpolant(xi,yi,sigma_n,'linear','nearest');
    ctr = (nodes(elems(:,1),:)+nodes(elems(:,2),:)+nodes(elems(:,3),:))/3;
    sigma_e = F(ctr(:,1),ctr(:,2));
    sigma_elem_all{i}=sigma_e;
    img = mk_image(fmdl,sigma_e);
    vh  = fwd_solve(img);
    allV{i}=vh.meas;
end
end
% ------------------------------------------------------------------------
function inv_model = buildEidorsDiffModel(n_elec)
% buildEidorsDiffModel  adjacent‑drive difference model (with background)
% 输入  n_elec  ——  电极数(同一环)
% 输出  inv_model —— 可直接用于 inv_solve(inv_model, v0, vt)

    % 1) 基础前向模型（粗网格即可加快重建）
    mdl = mk_common_model('c2c', n_elec);   % point‑elec, fast
    fmdl = mdl.fwd_model;

    % 2) 相邻激励 & 完整测量模式
    stim = mk_stim_patterns(n_elec, 1, [0,1], [0,1], {'do_redundant'}, 1);
    fmdl.stimulation = stim;

    % 3) 创建 inv_model
    inv_model                  = eidors_obj('inv_model', 'diff_adjacent');
    inv_model.fwd_model        = fmdl;
    inv_model.reconst_type     = 'difference';
    inv_model.solve            = @inv_solve_diff_pdipm;

    % 正则化与先验
    inv_model.hyperparameter.value = 1e-2;
    inv_model.RtR_prior            = @prior_laplace;
    inv_model.meas_constraint      = @meas_icov;

    % 4) **关键补丁** —— 提供背景图像
    img_bkg                     = mk_image(fmdl, 1);   % conductivity = 1
    inv_model.jacobian_bkgnd    = img_bkg;             % ★ 必须字段
end
% ------------------------------------------------------------------------
function img_list = doDiffReconstructionFromVoltage(allV, baseIdx, inv_model, t_vec)
% Visualize difference reconstruction for sanity‑check
T = numel(allV);  v0=allV{baseIdx};
figure('Name','Diff‑Recon','NumberTitle','off');
img_list=cell(T,1);
for k=1:T
    img = inv_solve(inv_model, v0, allV{k});
    img_list{k}=img;
    subplot(ceil(sqrt(T)),ceil(sqrt(T)),k);
    show_fem(img,[1,0]); axis equal tight off
    if nargin>3 && ~isempty(t_vec)
        title(sprintf('t=%.3f s', t_vec(k)));
    else
        title(sprintf('f%02d',k));
    end
end
end
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

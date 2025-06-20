 % ======================================================================
%  EIT_BuildDataset.m                                   2025‑06‑18 (rev)
% ----------------------------------------------------------------------
%  • 自动批量生成 EIT 触觉数据集（随机序列 + 常用固定组合）
%  • 全流程：COMSOL → Cut‑Plane Strain CSV → EIDORS → Electrode Voltage
%  • 额外输出段级标签 seg_info = {动作, 起始 t, 结束 t}
%  • 所有可调参数（位置 / 载荷幅值 / 作用范围 a / 运动时长 …）
%    统一在 cfg.PARAM_SWEEP.* 内声明，一键穷举
% ----------------------------------------------------------------------
%  依赖
%     – COMSOL 6.x LiveLink for MATLAB
%     – EIDORS v3.12‑ng
% ======================================================================
function EIT_BuildDataset()


%% ================= 0. 环境配置 ========================================
addpath('/mnt/disk2/software/comsol/mli');     
% ← COMSOL mli
run('~/projects/EIDORs/eidors-v3.12-ng/eidors/startup.m');          % EIDORS
import com.comsol.model.util.ModelUtil                              % COMSOL API

base_single = '/home/kairui/projects/EIT/Comsol6.2/EIT3D_SingleMotion.mph';
base_slide  = '/home/kairui/projects/EIT/Comsol6.2/EIT3D_Slide.mph';
out_root    = '/home/kairui/projects/EIT/Dataset';  if ~exist(out_root,'dir'), mkdir(out_root); end

%% ================ 1. 全局配置 ========================================
cfg = struct();

% —— 固定组合（最常用四连击） ————————————————
cfg.FIXED_PATTERNS = { {'Press','Slide:line','Tap','Punch'} };

% —— 随机序列参数 ————————————————————
cfg.N_RANDOM_SEQ = 8;         % 生成多少条随机序列
cfg.SEQ_LEN      = 4;         % 一条随机序列动作数
cfg.MAX_PER_KIND = 2;         % 同类动作在一条序列里出现上限
cfg.SINGLE_TYPES = {'Press','Tap','Punch'};

% —— 轨迹预设（Slide） ———————————————————
cfg.SLIDE_PRESET = { ...
    struct('type','line',  'k_x', 0.001,'k_y',  0.000), ...
    struct('type','line',  'k_x', 0.000,'k_y',  0.001), ...
    struct('type','circle','R',   0.003,'w', 2*pi/3) };

% —— 统一物理模型参数（不进 sweep） ——————————
cfg.MODEL_TAG = 'h2C';     % EIDORS 网格标签
cfg.GF       = 8;          % Gauge‑Factor
cfg.SIGMA0   = 1;

% —— 参数穷举 / 扫描表 ———————————————————
cfg.PARAM_SWEEP.x0  = linspace(-4e-3, 4e-3, 3);   % 位置 (m)
cfg.PARAM_SWEEP.y0  = linspace(-4e-3, 4e-3, 3);
cfg.PARAM_SWEEP.a   = [50, 100];                  % exp 衰减系数
cfg.PARAM_SWEEP.F   = [5e5, 1e6];                 % 力幅值 (Pa)
cfg.PARAM_SWEEP.Ts  = [3, 4, 5];                  % 单点动作时长 (s)
cfg.PARAM_SWEEP.Tv  = [2, 3, 4, 5];                  % Slide 动作时长 (s)

cfg.DT_single = 0.05;  cfg.DT_slide = 0.05;

rng('shuffle');

%% ================ 2. 构建序列规格列表 ================================
seq_specs = buildAllSpecs(cfg);   % cell 数组，每项一个规格
fprintf('[INFO] generate %d unsolved sequence pattern\n', numel(seq_specs));

%% ================ 3. 共享前向模型（EIDORS） ==========================
mdl_tmp     = mk_common_model(cfg.MODEL_TAG,16);   fmdl = mdl_tmp.fwd_model;

%% ================ 4. 循环求解每条序列 ===============================
for idx = 1:numel(seq_specs)
    spec = seq_specs{idx};
    fprintf('\n===== (%d/%d) %s =====\n', idx, numel(seq_specs), spec.id);

    [allV, t_vec, frame_tag, seg_info] = solveOneSpec(spec, cfg, base_single, base_slide, fmdl);

    save(fullfile(out_root, ['Seq_',spec.id,'.mat']), ...
        'allV','t_vec','frame_tag','seg_info','spec');
    fprintf('    ✔ save %s finish\n', spec.id);
end
end  % >>> main 结束

%% ====================================================================
%                       子函数区域                                      
% =====================================================================
function seq_specs = buildAllSpecs(cfg)
% • 固定组合先展开参数 sweep
% • 随机序列随后生成，参数同样走 sweep
seq_specs = {};

% —— 笛卡尔积生成参数组合 ——
paramGrid = allcomb(cfg.PARAM_SWEEP.x0, cfg.PARAM_SWEEP.y0, ...
                    cfg.PARAM_SWEEP.a,  cfg.PARAM_SWEEP.F);

%% —— 1) 固定序列 ——   
for p = 1:size(paramGrid,1)
    x0    = paramGrid(p,1);
    y0    = paramGrid(p,2);
    a_val = paramGrid(p,3);
    F_val = paramGrid(p,4);

    for fp = 1:numel(cfg.FIXED_PATTERNS)
        patt    = cfg.FIXED_PATTERNS{fp};
        motions = cellfun(@(tok) token2motion(tok,cfg), patt,'Uni',false);
        spec    = packSpec(x0,y0,a_val,F_val,motions,'FIX');
        seq_specs{end+1} = spec; %#ok<AGROW>
    end
end

% —— 2) 随机序列 ——   
for n = 1:cfg.N_RANDOM_SEQ
    % 随机抽一行参数
    ridx   = randi(size(paramGrid,1));
    x0     = paramGrid(ridx,1);
    y0     = paramGrid(ridx,2);
    a_val  = paramGrid(ridx,3);
    F_val  = paramGrid(ridx,4);

    motions = randomMotionList(cfg);
    spec    = packSpec(x0,y0,a_val,F_val,motions,'RND');
    seq_specs{end+1} = spec; %#ok<AGROW>
end

end

% ---------------------------------------------------------------------
function mlist = randomMotionList(cfg)
mlist = {}; cnt=0; cSlide=0; 
cSingle = containers.Map(cfg.SINGLE_TYPES, zeros(1, numel(cfg.SINGLE_TYPES)));

while cnt<cfg.SEQ_LEN
  if rand<.3 && cSlide<cfg.MAX_PER_KIND
      cSlide=cSlide+1; mv.kind='Slide'; mv.param=cfg.SLIDE_PRESET{randi(numel(cfg.SLIDE_PRESET))};
  else
      st = cfg.SINGLE_TYPES{randi(numel(cfg.SINGLE_TYPES))};
      if cSingle(st)<cfg.MAX_PER_KIND
          cSingle(st)=cSingle(st)+1; mv.kind=st; mv.param=[];
      else; continue; end
  end
  mlist{end+1}=mv; cnt=cnt+1; %#ok<AGROW>
end
end

% ---------------------------------------------------------------------
function motionspec = token2motion(tok,cfg)
if startsWith(tok,'Slide','IgnoreCase',true)
    parts = split(tok,':'); sub = lower(parts{min(2,end)});
    sw = struct('line',1,'sinex',2,'sine',2,'siney',3,'circle',4);
    motionspec.kind = 'Slide'; motionspec.param = cfg.SLIDE_PRESET{ sw.(sub) };
else
    motionspec.kind = tok; motionspec.param = [];
end
end

% ---------------------------------------------------------------------
function spec = packSpec(x0,y0,a_val,F_val,motions,prefix)
stamp = datestr(now,'yymmdd_HHMMSSFFF');
spec = struct('id',[prefix,'_',stamp,'_',char(java.util.UUID.randomUUID)], ...
              'x0',x0,'y0',y0,'a_val',a_val,'F',F_val,'motions',{motions});
end

%% ===================== 求解一条序列 =================================
function [allV,t_vec,frame_tag,seg_info] = solveOneSpec(spec,cfg,base_single,base_slide,fmdl)
allV={}; t_vec=[]; frame_tag={}; seg_info={};

com.comsol.model.util.ModelUtil.clear();

java.lang.System.gc();

for i=1:numel(spec.motions)
    mv = spec.motions{i};

   % — 时长与步长 —
    if strcmp(mv.kind,'Slide')
        T  = cfg.PARAM_SWEEP.Tv( randi(numel(cfg.PARAM_SWEEP.Tv)) );
        dt = cfg.DT_slide;
    else
        T  = cfg.PARAM_SWEEP.Ts( randi(numel(cfg.PARAM_SWEEP.Ts)) );
        dt = cfg.DT_single;
    end
    tlist   = sprintf('range(0,%g,%g)',dt,T);
    local_t = 0:dt:T;

      % — 载入 MPH —
    if strcmp(mv.kind,'Slide')
        model = mphload(base_slide,'tag','M');
    else
        model = mphload(base_single,'tag','M');
    end
    model.hist.disable();

    % —— 共用参数 ——
    model.param.set('x0', num2str(spec.x0));
    model.param.set('y0', num2str(spec.y0));
    model.param.set('a',  num2str(spec.a_val));
    model.param.set('force_scale', num2str(spec.F));

    % —— 动作具体设置 ——
    if strcmp(mv.kind,'Slide')
        cfgSlide(model,mv.param,spec.x0,spec.y0);
    else
        ftag=['F',lower(mv.kind)];
        expr=sprintf('- %s(t)*exp(-a*((x-x0)^2+(y-y0)^2))*%.0e',ftag,spec.F);
        model.physics('solid').feature('bndl1').set('FperArea',{expr,'0','0'});
    end
    model.study('std1').feature('time').set('tlist',tlist);
    model.study('std1').feature('time').set('useinitsol','off');
    model.study('std1').run;

    % —— 导出应变 ——
    csv=[tempname,'.csv']; ex=model.result.export.create('c','Data');
    ex.set('data','cpl1'); ex.set('expr',{'solid.ep1'}); ex.set('filename',csv); ex.run;
    strain=readCSV(csv); delete(csv);

    % —— 应变 → 电压 ——
    [~,Vcell]=strain2Voltage(strain,cfg.GF,cfg.SIGMA0,fmdl);

    % —— 汇总 ——
       if isempty(t_vec)
        t_offset = 0;
    else
        t_offset = t_vec(end);
    end

    t_vec    = [t_vec; local_t.'+t_offset];
    allV     = [allV;  Vcell];
    frame_tag= [frame_tag; repmat({mv.kind},numel(Vcell),1)];
    seg_info = [seg_info; {mv.kind, t_offset, t_offset+T}];

    com.comsol.model.util.ModelUtil.clear(); java.lang.System.gc();
end
end

%% ================= 工具函数 ==========================================
function cfgSlide(model,sp,x0,y0)
    switch sp.type
        case 'line'
            model.func('an1').set('expr',sprintf('%g+%g*t',x0,sp.k_x));
            model.func('an2').set('expr',sprintf('%g+%g*t',y0,sp.k_y));
        case 'circle'
            model.func('an1').set('expr',sprintf('%g+%g*cos(%g*t)',x0,sp.R,sp.w));
            model.func('an2').set('expr',sprintf('%g+%g*sin(%g*t)',y0,sp.R,sp.w));
        case 'sine'
            if sp.axis=='x'
               model.func('an1').set('expr',sprintf('%g+%g*sin(%g*t)',x0,sp.A,sp.w));
               model.func('an2').set('expr',sprintf('%g',y0));
            else
               model.func('an1').set('expr',sprintf('%g',x0));
               model.func('an2').set('expr',sprintf('%g+%g*sin(%g*t)',y0,sp.A,sp.w));
            end
    end
end

function data = readCSV(f)
    fid=fopen(f); for i=1:4,fgetl(fid); end
    N=sscanf(fgetl(fid),'%% Nodes,%d'); T=sscanf(fgetl(fid),'%% Expressions,%d');
    for i=1:3,fgetl(fid); end; raw=zeros(N,3+T);
    for n=1:N, raw(n,:)=sscanf(fgetl(fid),'%f,').'; end; fclose(fid);
    x=raw(:,1); y=raw(:,2); vs=raw(:,4:end);
    data=zeros(T,N,3);
    for t=1:T, data(t,:,1)=x.'; data(t,:,2)=y.'; data(t,:,3)=vs(:,t).'; end
end



%%=======================================================================%%
% allcomb  ——  n 维笛卡儿积（每行是一组组合）
%    C = allcomb(A,B,...) , 其中 A,B... 均为向量或标量
%    例： allcomb([1 2],[3 4])  →  [1 3; 1 4; 2 3; 2 4]
%%=======================================================================%%
function C = allcomb(varargin)
    assert(nargin >= 1, 'allcomb: 至少需要一个输入');
    % 1) 把所有输入拉成列向量
    grids = cellfun(@(v)v(:), varargin, 'uni', false);
    % 2) ndgrid 展开
    [grids{:}] = ndgrid(grids{:});
    % 3) 拼到一起，reshape 成 N×k
    k = numel(grids);
    C = reshape(cat(k+1, grids{:}), [], k);
end

%====================================================================
%  全新版本 —— 先去重，再插值
%====================================================================
function [sigma_elem_all, Vcell] = strain2Voltage(bigArr, GF, sigma0, fmdl)

    % ---------- 预计算三角形重心 ----------
    persistent ctr
    if isempty(ctr)
        nodes = fmdl.nodes; elems = fmdl.elems;
        ctr = ( nodes(elems(:,1),:) + nodes(elems(:,2),:) + nodes(elems(:,3),:) ) / 3;
    end

    [T,~,~] = size(bigArr);
    sigma_elem_all = cell(T,1);  Vcell = cell(T,1);

    for t = 1:T
        % === ① 取出该帧节点坐标 / 应变 ===
        xi  = bigArr(t,:,1)';           % N ×1
        yi  = bigArr(t,:,2)';           % N ×1
        eps = bigArr(t,:,3)';           % N ×1

        % === ② “坐标 → 行向量” 形成 key，unique 去重 ===
        XY   = [xi yi];                 % N ×2
        [XYu,~,idx] = unique(XY,'rows');% M ×2,  idx:N→M
        epsu = accumarray(idx, eps, [], @mean);   % 对重复点求平均

        % === ③ 电导率节点值 ===
        sigma_n = sigma0 .* (1 - GF .* epsu);

        % === ④ 插值到单元中心 → 元电导率 ===
        F = scatteredInterpolant(XYu(:,1), XYu(:,2), sigma_n, ...
                                 'linear', 'nearest');
        sigma_e = F(ctr(:,1), ctr(:,2));           % (#elems) ×1

        % === ⑤ EIDORS 前向求解 ===
        img   = mk_image(fmdl, sigma_e);
        vh    = fwd_solve(img);

        % === ⑥ 保存 ===
        sigma_elem_all{t} = sigma_e;
        Vcell{t}          = vh.meas;
    end
end

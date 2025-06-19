% ======================================================================
%  EIT_BuildDataset.m                                      2025‑06‑17
% ----------------------------------------------------------------------
%  • One‑click dataset builder for Electrical‑Impedance‑Tactile project
%  • Generates *both* random and user‑preset motion sequences in batch
%  • Each sequence ⇒ COMSOL → strain CSV ⇒ EIDORS → electrode voltage
%  • All parameters are configured **once, below**, no runtime inputs
% ----------------------------------------------------------------------
%  Requirements
%     – COMSOL LiveLink for MATLAB                   (tested on 6.2)
%     – EIDORS toolbox                               (tested on 3.12‑ng)
% ======================================================================

function EIT_BuildDataset()
%% ================= 0. Environment =====================================
addpath('/mnt/disk2/software/comsol/mli');                    % ← edit if needed
eidors_startup = '/home/kairui/projects/EIDORs/eidors-v3.12-ng/eidors/startup.m';
run(eidors_startup);                                          % EIDORS
import com.comsol.model.util.ModelUtil                        % COMSOL API

% -------- basic paths --------
base_single = '/home/kairui/projects/EIT/Comsol6.2/EIT3D_SingleMotion.mph';
base_slide  = '/home/kairui/projects/EIT/Comsol6.2/EIT3D_Slide.mph';
out_root    = '/home/kairui/projects/EIT/Dataset';
if ~exist(out_root,'dir'); mkdir(out_root); end

%% ================ 1. Global configuration struct ======================
cfg = struct();

% ---- sequence pool size ----
cfg.N_RANDOM_SEQ   = 2;        % how many *random* sequences to create
cfg.FIXED_SEQS     = { ...     % cell array of *fixed* motion token lists
    {'Press','Slide:line','Tap','Punch'}, ...
    {'Slide:circle','Press','Slide:sineX','Tap'} };

% ---- motion composition ----
cfg.SEQ_LEN        = 4;        % length for *random* sequences
cfg.MAX_PER_KIND   = 2;        % per‑kind upper bound inside a sequence
cfg.SINGLE_TYPES   = {'Press','Tap','Punch'};

% ---- slide trajectory presets ----
cfg.SLIDE_PARAM_SET = { ...
    struct('type','line',  'k_x', 0.001,'k_y',  0.000), ...
    struct('type','line',  'k_x', 0.001,'k_y', -0.001), ...
    struct('type','line',  'k_x', 0.000,'k_y',  0.001), ...
    struct('type','sine',  'axis','x','A', 0.003,'w', 2*pi/3), ...
    struct('type','sine',  'axis','y','A', 0.003,'w', 2*pi/3), ...
    struct('type','circle','R',   0.003,'w', 2*pi/3)};

% ---- mechanical / physical constants ----
cfg.x0_range  = [-0.005, 0.005];  % m
cfg.y0_range  = [-0.005, 0.005];
cfg.a_val     = 100;              % exp(‑a·r²)
cfg.force_val = 1e6;              % N → Pa scaling

% ---- simulation time control (randomised) ----
cfg.TLEN_RANGE_single = [2, 4];   % s  (Press/Tap/Punch)
cfg.TLEN_RANGE_slide  = [2, 4];   % s  (Slide)
cfg.DT_single         = 0.05;     % s  step size for single
cfg.DT_slide          = 0.05;     % s  step size for slide

% ---- strain → voltage ----
cfg.GF        = 8;
cfg.SIGMA0    = 1;
cfg.MODEL_TAG = 'h2C';            % 4096‑elem, 16‑elec CEM 2‑D model

rng('shuffle');                   % global random seed

%% ================ 2. Build specification lists ========================
spec_rand  = buildRandomSpec(cfg);
spec_fixed = buildFixedSpec(cfg);
all_spec   = [spec_rand; spec_fixed];     % ← vertical concatenation to
                                          %    avoid size mismatch
fprintf('[INFO] Total sequences – random: %d  fixed: %d\n', ...
        numel(spec_rand), numel(spec_fixed));

%% ================ 3. Shared EIDORS forward model ======================
mdl_tmp      = mk_common_model(cfg.MODEL_TAG,16);
fmdl_global  = mdl_tmp.fwd_model;

%% ================ 4. Generate each sequence ===========================
for k = 1:numel(all_spec)
    seq = all_spec{k};
    fprintf('\n===== Generating sequence %d / %d =====\n', k, numel(all_spec));
    seq_id = datestr(now,'yyyymmdd_HHMMSSFFF');   % unique id

    [allV, t_vec, motion_tag] = generateSequence(seq, cfg, ...
                                   base_single, base_slide, fmdl_global);

    save(fullfile(out_root, ['Seq_',seq_id,'.mat']), ...
         'allV','t_vec','motion_tag','seq');
    fprintf('[OK] Saved MAT → %s\n', seq_id);
end

end  % ================= END MAIN =======================================

%% ======================================================================
%                         SUB‑FUNCTIONS                                    
% ======================================================================
function specList = buildRandomSpec(cfg)
% buildRandomSpec  → cell{N_RANDOM_SEQ} each contains .x0 .y0 .motions
specList = cell(cfg.N_RANDOM_SEQ,1);
for n = 1:cfg.N_RANDOM_SEQ
    x0 = randRange(cfg.x0_range);
    y0 = randRange(cfg.y0_range);

    motions = {};  cntSlide = 0;
    cntSingle = containers.Map(cfg.SINGLE_TYPES, zeros(1,numel(cfg.SINGLE_TYPES)));
    while numel(motions) < cfg.SEQ_LEN
        if rand < .28 && cntSlide < cfg.MAX_PER_KIND   % choose Slide 28‑%
            cntSlide = cntSlide + 1;
            mv.kind  = 'Slide';
            mv.param = cfg.SLIDE_PARAM_SET{ randi(numel(cfg.SLIDE_PARAM_SET)) };
            motions{end+1} = mv; %#ok<AGROW>
        else
            stype = cfg.SINGLE_TYPES{ randi(numel(cfg.SINGLE_TYPES)) };
            if cntSingle(stype) < cfg.MAX_PER_KIND
                cntSingle(stype) = cntSingle(stype) + 1;
                mv.kind  = stype;  mv.param = [];
                motions{end+1} = mv; %#ok<AGROW>
            end
        end
    end
    specList{n} = struct('x0',x0,'y0',y0,'motions',{motions});
end
end

% ----------------------------------------------------------------------
function specList = buildFixedSpec(cfg)
% buildFixedSpec – turn token lists into spec structs
specList = cell(numel(cfg.FIXED_SEQS),1);
for i = 1:numel(cfg.FIXED_SEQS)
    tokens = cfg.FIXED_SEQS{i};
    x0 = randRange(cfg.x0_range);
    y0 = randRange(cfg.y0_range);
    motions = cellfun(@(s)parseOneToken(s,cfg), tokens, 'Uni',false);
    specList{i} = struct('x0',x0,'y0',y0,'motions',{motions});
end
end

function mv = parseOneToken(tok,cfg)
if startsWith(tok,'Slide','IgnoreCase',true)
    mv.kind = 'Slide';
    parts = split(tok,':');
    subtype = lower(parts{min(end,2)});
    switch subtype
        case {'line'};   mv.param = cfg.SLIDE_PARAM_SET{1};
        case {'sinex','sine'}; mv.param = cfg.SLIDE_PARAM_SET{4};
        case {'siney'};  mv.param = cfg.SLIDE_PARAM_SET{5};
        case {'circle'}; mv.param = cfg.SLIDE_PARAM_SET{6};
        otherwise; error('Unknown Slide subtype %s',subtype);
    end
else
    mv.kind  = strtrim(tok);
    mv.param = [];
end
end

% ----------------------------------------------------------------------
function [allV, t_vec, motion_tag] = generateSequence(seq, cfg, ...
                        base_single, base_slide, fmdl_global)
allV = {};  t_vec = [];  motion_tag = {};

for mIdx = 1:numel(seq.motions)
    mv = seq.motions{mIdx};
    % ===== time list (random length) =====
    if strcmp(mv.kind,'Slide')
        T_len = randRange(cfg.TLEN_RANGE_slide); dt = cfg.DT_slide;
    else
        T_len = randRange(cfg.TLEN_RANGE_single); dt = cfg.DT_single;
    end
    tlistStr = sprintf('range(0,%g,%g)', dt, T_len);

    % ===== load MPH =====
    if strcmp(mv.kind,'Slide'); model = mphload(base_slide);model.label('WorkingModel'); else; model = mphload(base_single); end

    % ===== common parameters =====
    model.param.set('x0', num2str(seq.x0));
    model.param.set('y0', num2str(seq.y0));
    model.param.set('a',  num2str(cfg.a_val));
    model.param.set('force_scale', num2str(cfg.force_val));

    % ===== motion‑specific setup =====
    if strcmp(mv.kind,'Slide')
        setSlideTrajectory(model, mv.param, seq.x0, seq.y0);
    else
        ftag = ['F', lower(mv.kind)];
        expr = sprintf('- %s(t)*exp(-a*((x-x0)^2+(y-y0)^2))*%.0e', ...
                       ftag, cfg.force_val);
        model.physics('solid').feature('bndl1').set('FperArea',{expr,'0','0'});
    end
    model.study('std1').feature('time').set('tlist', tlistStr);
    model.study('std1').feature('time').set('useinitsol','off');
    model.study('std1').run;

    % ===== export strain CSV =====
    tmp_csv = [tempname,'.csv'];
    ex = model.result.export.create('cStrain','Data');
    ex.set('data','cpl1'); ex.set('expr',{'solid.ep1'}); ex.set('filename',tmp_csv); ex.run;
    strain = readCsv(tmp_csv); delete(tmp_csv);

    % ===== strain → voltage =====
    [~, Vcell] = strain2Voltage(strain, cfg.GF, cfg.SIGMA0, fmdl_global);

    % ===== concatenate =====
    local_t = evalT(tlistStr);
   if isempty(t_vec)
    offset = 0;
   else
    offset = t_vec(end);
   end
    t_vec       = [t_vec; local_t + offset]; %#ok<AGROW>
    allV        = [allV;  Vcell];            %#ok<AGROW>
    motion_tag  = [motion_tag; repmat({mv.kind}, numel(Vcell),1)]; %#ok<AGROW>

    com.comsol.model.util.ModelUtil.remove('WorkingModel');
    java.lang.System.gc(); 
    pause(0.05);
end
end

%% ======================================================================
%                           TOOL FUNCTIONS                                
% ======================================================================
function v = randRange(r); v = r(1) + (r(2)-r(1))*rand; end
function t = evalT(str); n=sscanf(str,'range(%f,%f,%f)'); t=(n(1):n(2):n(3)).'; end

function setSlideTrajectory(model, sp, x0, y0)
    switch sp.type
        case 'line'
            model.func('an1').set('expr', sprintf('%g+%g*t', x0, sp.k_x));
            model.func('an2').set('expr', sprintf('%g+%g*t', y0, sp.k_y));
        case 'sine'
            if sp.axis=='x'
                model.func('an1').set('expr', sprintf('%g+%g*sin(%g*t)', x0, sp.A, sp.w));
                model.func('an2').set('expr', sprintf('%g', y0));
            else
                model.func('an1').set('expr', sprintf('%g', x0));
                model.func('an2').set('expr', sprintf('%g+%g*sin(%g*t)', y0, sp.A, sp.w));
            end
        case 'circle'
            model.func('an1').set('expr', sprintf('%g+%g*cos(%g*t)', x0, sp.R, sp.w));
            model.func('an2').set('expr', sprintf('%g+%g*sin(%g*t)', y0, sp.R, sp.w));
        otherwise
            error('Unknown slide type %s', sp.type);
    end
end

function strain = readCsv(csvFile)
    fid=fopen(csvFile,'r'); for i=1:4, fgetl(fid); end
    N=sscanf(fgetl(fid),'%% Nodes,%d'); T=sscanf(fgetl(fid),'%% Expressions,%d');
    for i=1:3, fgetl(fid); end
    raw=zeros(N,3+T);
    for n=1:N, raw(n,:)=sscanf(fgetl(fid),'%f,').'; end; fclose(fid);
    x=raw(:,1); y=raw(:,2); vs=raw(:,4:end);
    strain=zeros(T,N,3);
    for t=1:T
        strain(t,:,1)=x.'; strain(t,:,2)=y.'; strain(t,:,3)=vs(:,t).';
    end
end

function [sigma_elem_all, Vcell] = strain2Voltage(bigArr, GF, sigma0, fmdl)
    [T,~,~]=size(bigArr); nodes=fmdl.nodes; elems=fmdl.elems;
    sigma_elem_all=cell(T,1); Vcell=cell(T,1);
    for i=1:T
        xi=bigArr(i,:,1)'; yi=bigArr(i,:,2)'; eps=bigArr(i,:,3)';
        sigma_n = sigma0 .* (1 - GF .* eps);
        F = scatteredInterpolant(xi,yi,sigma_n,'linear','nearest');
        ctr = (nodes(elems(:,1),:)+nodes(elems(:,2),:)+nodes(elems(:,3),:))/3;
        sigma_e = F(ctr(:,1),ctr(:,2));
        sigma_elem_all{i}=sigma_e;
        img = mk_image(fmdl,sigma_e); vh = fwd_solve(img);
        Vcell{i} = vh.meas;
    end
end

function r = tern(cond,x,y); if cond; r=x; else; r=y; end; end

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

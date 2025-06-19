%% =====  Environment  =====
addpath('/mnt/disk2/software/comsol/mli');      % <- change if needed
import com.comsol.model.util.ModelUtil          % LiveLink

%% =====  User Config  =====
base_single = '/home/kairui/projects/EIT/Comsol6.2/EIT3D_SingleMotion.mph';           % mph for press/tap/punch
base_slide  = '/home/kairui/projects/EIT/Comsol6.2/EIT3D_Slide.mph';                  % mph for slide
out_dir     = '/home/kairui/projects/EIT/RandCompositeOut';
if ~exist(out_dir,'dir'), mkdir(out_dir); end

% ----   random sequence settings   ----
SEQ_LEN        = 4;                             % how many motions in one run
SINGLE_TYPES   = {'Press','Tap','Punch'};       % tags used in mph
SLIDE_PARAM_SET = {                           	% candidate slide param struct set
    struct('k_x',0.001,'k_y', 0.000), ...
    struct('k_x',0.001,'k_y',-0.001), ...
    struct('k_x',0.000,'k_y', 0.001) };

% ----   constant parameters   ----
x0_range   = [-0.005, 0.005];                   % mm → m (modify to your device)
y0_range   = [-0.005, 0.005];
a_val      = 100;                               % decay radius
force_val  = 1e9;                               % same force for every motion
tlist_single = 'range(0,0.03,3)';
tlist_slide  = 'range(0,0.03,3)';

rng('shuffle');                                 % random seed

%% =====  Generate one composite  =====
% 1) randomly pick base (x0,y0)
x0 = randRange(x0_range);
y0 = randRange(y0_range);

% 2) randomly build a motion list, first / last 均可能是 slide 或 single
motions = cell(SEQ_LEN,1);
for i = 1:SEQ_LEN
    if rand < .4                      % ~40% chance choose slide
        motions{i}.kind = 'Slide';
        % draw random slide param from list
        sp = SLIDE_PARAM_SET{ randi(numel(SLIDE_PARAM_SET)) };
        motions{i}.param = sp;
    else                              % single
        stype = SINGLE_TYPES{ randi(numel(SINGLE_TYPES)) };
        motions{i}.kind  = stype;
        motions{i}.param = [];        %#ok<AGROW>
    end
end
fprintf("◆ Sequence:\n");
cellfun(@(m)fprintf("   - %s\n",m.kind),motions);

%% =====  Loop and solve  =====
case_id  = datestr(now,'yyyymmdd_HHMMSS');
all_strain = [];        % concatenated [T,N,3]
all_time   = [];        % global time vector
motion_tag = {};        % store label for each frame

for m = motions'
    mk = m{1}.kind;

    switch mk
        case 'Slide'
            mph_file = base_slide;
        otherwise   % press/tap/punch
            mph_file = base_single;
    end
    model = mphload(mph_file);

    % ---------- set shared parameters ----------
    model.param.set('x0', num2str(x0));
    model.param.set('y0', num2str(y0));
    model.param.set('a',  num2str(a_val));
    model.param.set('force_scale', num2str(force_val));

    % ---------- specific for motion ----------
    if strcmp(mk,'Slide')
        % slide uses analytic functions an1, an2: modify slope
        sp = m{1}.param;
        model.func('an1').set('expr', sprintf('%g+%g*t',x0,sp.k_x));
        model.func('an2').set('expr', sprintf('%g+%g*t',y0,sp.k_y));
        model.study('std1').feature('time').set('tlist', tlist_slide);
    else
        % choose boundary load function tag automatically
        ftag = ['F', lower(mk)];      % e.g. 'Press' -> 'Fpress'
        expr = sprintf('- %s(t)*exp(-a*((x-x0)^2+(y-y0)^2))*%.0e', ...
                       ftag, force_val);
        model.physics('solid').feature('bndl1') ...
             .set('FperArea',{expr,'0','0'});
        model.study('std1').feature('time').set('tlist', tlist_single);
    end
    model.study('std1').feature('time').set('useinitsol','off');

    % ---------- solve ----------
    model.study('std1').run;
    fprintf('   ✔ %s solved\n',mk);

    % ---------- export strain csv ----------
    tmp_csv = fullfile(out_dir,'tmp.csv');
    if any(string(model.result.export.tags)=="expCsv")
        model.result.export.remove('expCsv');
    end
    ex = model.result.export.create('expCsv','Data');
    ex.set('data','cpl1');
    ex.set('expr',{'solid.ep1'});
    ex.set('filename',tmp_csv); ex.run;

    strain = readCsv(tmp_csv); delete(tmp_csv);
    stepN  = size(strain,1);                 % #time frames of this motion
    all_strain = cat(1, all_strain, strain); % concat along time dim

    % time vector for this motion
    if strcmp(mk,'Slide')
        tvec = evalTlist(tlist_slide);
    else
        tvec = evalTlist(tlist_single);
    end
    % ----------- safe append time vector -----------
if isempty(all_time)
    offset = 0;
else
    offset = all_time(end);   % last time stamp as offset
end
all_time = [all_time;  tvec(:) + offset];   % concatenate

    % tag
    motion_tag = [motion_tag; repmat({mk},stepN,1)]; %#ok<AGROW>

    ModelUtil.remove('model');
    java.lang.System.gc(); pause(0.1);
end

%% =====  Save mat =====
save_file = fullfile(out_dir, sprintf('Composite_%s.mat',case_id));
save(save_file,'all_strain','all_time','motion_tag','x0','y0', ...
               'a_val','force_val','motions');
fprintf("\n★ Saved composite MAT → %s\n",save_file);

%% =====  Helper functions =====
function v = randRange(r)
% uniform random value in [min,max]
v = r(1) + (r(2)-r(1))*rand;
end

function t = evalTlist(str)
% convert 'range(start,step,end)' to numeric vector
num = sscanf(str,'range(%f,%f,%f)');
t   = num(1):num(2):num(3);
t   = t(:);
end
%% ------------------------------------------------------------------------
function strain = readCsv(csvFile)
% readCsv  Read Cut‑Plane CSV exported from COMSOL
%   Return 3‑D array: [T, N, 3]  (T = time frames, N = nodes)
%       dim3 = {x, y,   strain}
%
%   CSV layout (COMSOL default):
%   row1‑4 : header
%   row5   : "% Nodes, <N>"
%   row6   : "% Expressions, <T>"
%   row7‑8 : header
%   row9   : "x (m), y (m), z (m), solid.ep1 (1), solid.ep1 (2) …"
%   next N rows: data

    fid = fopen(csvFile,'r');
    if fid<0
        error("Cannot open %s",csvFile);
    end

    for i = 1:4, fgetl(fid); end           % skip 4 header lines
    N = sscanf(fgetl(fid),'%% Nodes,%d');   % #nodes
    T = sscanf(fgetl(fid),'%% Expressions,%d');  % #time frames
    fgetl(fid); fgetl(fid); fgetl(fid);     % skip 3 more lines (titles)

    raw = zeros(N,3+T);                     % x,y,z, val1..valT
    for n = 1:N
        raw(n,:) = sscanf(fgetl(fid),'%f,').';
    end
    fclose(fid);

    x  = raw(:,1);
    y  = raw(:,2);
    vs = raw(:,4:end);                      % [N × T]

    strain = zeros(T,N,3);
    for t = 1:T
        strain(t,:,1) = x.';                % x coord
        strain(t,:,2) = y.';                % y coord
        strain(t,:,3) = vs(:,t).';          % strain value
    end
end

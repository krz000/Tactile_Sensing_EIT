%% ====   环境&库   ====
addpath('/mnt/disk2/software/comsol/mli');
import com.comsol.model.util.ModelUtil

%% ====   用户参数   ====
base_model = '/home/kairui/projects/EIT/Comsol6.2/EIT3D_Slide.mph';
out_dir    = '/home/kairui/projects/EIT/SlideTraj';
if ~exist(out_dir,'dir'), mkdir(out_dir); end

tlist       = 'range(0,0.03,3)';
a_val       = 100;        % exp 衰减
force_scale = 1e9;        % 幅值

trajDefs = {                       % 轨迹名 | 参数结构 | x(t) | y(t)
 'Line', struct('x0',[0 2e-3],'y0',0,'k_x',[5e-4 1e-3],'k_y',[5e-4 -1e-3]), ...
          @(p)sprintf('%g+%g*t',p.x0,p.k_x), ...
          @(p)sprintf('%g+%g*t',p.y0,p.k_y);

 'Circle', struct('x0',0,'y0',0,'R',[3e-3 4e-3],'w',[15 30]), ...
          @(p)sprintf('%g+%g*cos(%g*t)',p.x0,p.R,p.w), ...
          @(p)sprintf('%g+%g*sin(%g*t)',p.y0,p.R,p.w);

 'SineX', struct('x0',0,'y0',[0 2e-3],'A',[2e-3 3e-3],'w',[20 35]), ...
          @(p)sprintf('%g+%g*sin(%g*t)',p.x0,p.A,p.w), ...
          @(p)sprintf('%g',p.y0);
};

%% ====   主循环   ====
cid = 0;
for it = 1:size(trajDefs,1)
    name  = trajDefs{it,1};
    gS    = trajDefs{it,2};
    fxFun = trajDefs{it,3};
    fyFun = trajDefs{it,4};

    % --------- 自动展开全部参数组合 ----------
    fld  = fieldnames(gS);
    vecs = cellfun(@(f) gS.(f)(:), fld,'uni',0);
    [vecs{:}] = ndgrid(vecs{:});               % 自适应字段数
    M   = numel(vecs{1});
    flat = cellfun(@(v) reshape(v,[],1), vecs,'uni',0);
    mat  = [flat{:}];                          % M × F
    paramList = arrayfun(@(i) ...
               cell2struct(num2cell(mat(i,:)), fld, 2), ...
               1:M, 'uni',1).';

    % ------------- 对每组参数求解 -------------
    for p = paramList.'
        cid = cid+1;
        fprintf('\n[%03d] %s  %s\n',cid,name,jsonencode(p));
       

        mdl = mphload(base_model);

        mdl.func('an1').set('expr', fxFun(p));
        mdl.func('an2').set('expr', fyFun(p));
        mdl.param.set('a', num2str(a_val));
        mdl.param.set('force_scale',num2str(force_scale));

        mdl.study('std1').feature('time').set('tlist',tlist);
        mdl.study('std1').feature('time').set('useinitsol','off');
        mdl.study('std1').run;

        tmp = sprintf('%s/%s_%03d_tmp.csv',out_dir,name,cid);
        if any(string(mdl.result.export.tags)=="expCsv")
            mdl.result.export.remove('expCsv');
        end
        ex = mdl.result.export.create('expCsv','Data');
        ex.set('data','cpl1');
        ex.set('expr',{'solid.ep1'});
        ex.set('filename',tmp);
        ex.run;

        strain = readCsv(tmp); delete(tmp);
        save(fullfile(out_dir,sprintf('%s_%03d.mat',name,cid)), ...
             'strain','name','p','a_val','force_scale','tlist');

        fprintf('       ✔ Saved (%d×%d)\n',size(strain,2),size(strain,1));
        ModelUtil.remove('model');
        java.lang.System.gc(); pause(0.2);
    end
end
fprintf('\n🎉 Done! All files in %s\n',out_dir);

%% ====== CSV 读取 ======
function A = readCsv(f)
    fid=fopen(f); for i=1:4,fgetl(fid);end
    N=sscanf(fgetl(fid),'%% Nodes,%d'); T=sscanf(fgetl(fid),'%% Expressions,%d');
    fgetl(fid);fgetl(fid);fgetl(fid);
    D=zeros(N,3+T);
    for n=1:N, D(n,:)=sscanf(fgetl(fid),'%f,').'; end, fclose(fid);
    x=D(:,1);y=D(:,2);v=D(:,4:end); A=zeros(T,N,3);
    for t=1:T, A(t,:,1)=x.'; A(t,:,2)=y.'; A(t,:,3)=v(:,t).'; end
end

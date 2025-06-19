%% =====  EIT 复合动作数据生成器（多参数 + 复合动作） =====
import com.comsol.model.*
import com.comsol.model.util.*

base_model = '/home/kairui/projects/EIT/Comsol6.2/EIT3D_SingleMotion.mph';
actions    = {'Press','Fpress'; 'Tap','Ftap';'Punch','Fpunch'};

% —— 在这里定义参数列表 ——  
x0_list         = [0.00, 0.01, 0.02];
y0_list         = [0.00, 0.01, 0.02];
a_list          = [ 100,  200,  300];
force_scale_list= [1e9, 2e9];

tlist   = 'range(0,0.03,3)';
outDir  = '/home/kairui/projects/EIT/SingleMotionCompositeActions';
if ~exist(outDir,'dir'), mkdir(outDir); end

% ---------- 加载模型 & 锁定 tag ----------
model   = mphload(base_model);
studyArr= string(model.study().tags());    stdTag  = char(studyArr(1));
stepArr = string(model.study(stdTag).feature.tags());  timeTag = char(stepArr(1));
model.study(stdTag).feature(timeTag).set('tlist',tlist);

dsetArr = string(model.result.dataset.tags());
cpTag   = char(dsetArr(contains(lower(dsetArr),'cpl')));

fprintf('>>> Study:%s  Step:%s  CutPlane:%s\n',stdTag,timeTag,cpTag);

% ---------- 预置导出节点 ----------
if any(string(model.result.export.tags())=="expCsv")
    model.result.export.remove('expCsv');
end
ex = model.result.export.create('expCsv','Data');
ex.set('data',cpTag); ex.set('expr',{'solid.ep1'});

% ---------- 参数循环 + 主循环 ----------
for ix = 1:numel(x0_list)
  for iy = 1:numel(y0_list)
    for ia = 1:numel(a_list)
      for ifs = 1:numel(force_scale_list)

        % —— 当前参数 ——  
        x0          = x0_list(ix);
        y0          = y0_list(iy);
        a           = a_list(ia);
        force_scale = force_scale_list(ifs);

        fprintf('\n>> Parameter Combination: x0=%.3f y0=%.3f a=%d f=%.1e\n',...
                x0,y0,a,force_scale);

        % 清空合并容器
        compData  = [];
        compLabel = "";

        % —— 复合动作循环 ——  
        for k = 1:size(actions,1)
          act  = actions{k,1};
          ftag = actions{k,2};
          fprintf('  [%d/%d] Motion %-5s\n',k,size(actions,1),act);

          % 1) 更新参数 & 载荷
          model.param.set('x0',num2str(x0));
          model.param.set('y0',num2str(y0));
          model.param.set('a', num2str(a));
          expr = sprintf('- %s(t)*exp(-a*((x-x0)^2+(y-y0)^2))*%.0e',...
                         ftag,force_scale);
          model.physics('solid').feature('bndl1')...
               .set('FperArea',{expr,'0','0'});

          % 2) 强制零初值
          model.study(stdTag).feature(timeTag)...
               .set('useinitsol','off');

          % 3) 运行
          model.study(stdTag).run();

          % 4) 导出 → 读 CSV → 合并
          tmp = fullfile(outDir,sprintf('%s_x%.3f_y%.3f_a%d_f%.1e_tmp.csv',...
                                        act,x0,y0,a,force_scale));
          ex.set('filename',tmp); ex.run();
          arr = readCsv(tmp); delete(tmp);

          if k==1
            compData  = arr;
            compLabel = sprintf('%s',act);
          else
            compData  = cat(1,compData,arr);
            compLabel = [compLabel '+' act];
          end
          fprintf('    ✔ %s Composite Done ( %d Steps)\n',...
                  act, size(compData,1));
        end

        % —— 保存本参数组合的结果 ——  
        outFile = fullfile(outDir,...
                   sprintf('Composite_%s_x%.3f_y%.3f_a%d_f%.1e.mat',...
                           compLabel,x0,y0,a,force_scale));
        save(outFile,'compData','x0','y0','a','force_scale','compLabel');
        fprintf('  [√] saved to %s\n', outFile);

      end
    end
  end
end

%% ---------- 辅助函数 ----------
function A = readCsv(f)
  fid = fopen(f); for i=1:4, fgetl(fid); end
  N   = sscanf(fgetl(fid),'%% Nodes,%d');
  t   = sscanf(fgetl(fid),'%% Expressions,%d');
  fgetl(fid); fgetl(fid); fgetl(fid);
  D   = zeros(N,3+t);
  for n=1:N, D(n,:) = sscanf(fgetl(fid),'%f,').'; end
  fclose(fid);
  x = D(:,1);  y = D(:,2);  v = D(:,4:end);
  A = zeros(t,N,3);
  for i=1:t
    A(i,:,1) = x.';
    A(i,:,2) = y.';
    A(i,:,3) = v(:,i).';
  end
end

function EIT_ParamFit_Strain()
% Batch-fit Young’s modulus (E) & Poisson’s ratio (nu) by matching ΔR/R0.
% Verified on COMSOL 6.2 LiveLink for MATLAB.
%
% Steps
%   1.  Define paths / constants (中文：路径与常量)
%   2.  Double-loop over E × nu
%   3.  Per curve: update material ↔ force, run study, export strain
%   4.  Convert ε → ΔR/R0, compute RMSE, save plot
%   5.  Append summary.csv
%
% -------------------------------------------------------------------------

%% 1 ────────── PATHS & GLOBAL CONSTANTS (全局常量) ─────────────────────
baseModel = '/home/kairui/projects/EIT/Comsol6.2/EIT3D_SingleMotion.mph';

forceTxt  = { ...
    '/home/kairui/projects/EIT/Force_Curve/deformation/force_light.txt', ...
    '/home/kairui/projects/EIT/Force_Curve/deformation/force_strong.txt', ...
    '/home/kairui/projects/EIT/Force_Curve/deformation/force_longholdStrong.txt'};

strainTxt = { ...
    '/home/kairui/projects/EIT/Force_Curve/deformation/deformation_light.txt', ...
    '/home/kairui/projects/EIT/Force_Curve/deformation/deformation_strong.txt', ...
    '/home/kairui/projects/EIT/Force_Curve/deformation/deformation_longholdStrong.txt'};

x0 = 0.00;  y0 = 0.00;         % 提取中心
R  = 0.003;                    % 平均半径
aVal       = 100;              % 载荷分布 exp(-a·r²)
forceScale = 1e3;              % 力 → 压强缩放
GF         = 8;                % gauge factor
tlistStr   = 'range(0,0.05,4.2)';   % 固定求解时间

EList  = [160e9, 200e9, 240e9, 280e9, 360e9];
nuList = [0.22, 0.25, 0.275, 0.30];

outDir = fullfile(pwd,'fit_results');
if ~exist(outDir,'dir'), mkdir(outDir); end

import com.comsol.model.*
import com.comsol.model.util.*

%% 2 ────────── PARAMETRIC LOOP (双重循环) ─────────────────────────────
for E = EList
    for nu = nuList
        fprintf('\n===  E = %.2e  |  ν = %.3f  ===\n',E,nu);
        rmseAll = NaN(1,numel(forceTxt));   % 记录各曲线 RMSE

        %% 针对每条 (force, strain) 曲线
        for k = 1:numel(forceTxt)
            fprintf('  » Case #%d  ...  ',k);

            % ---------- 参照 ΔR/R0 ----------
            refMat = readmatrix(strainTxt{k});
            if size(refMat,2)<2 || all(~isfinite(refMat(:,2)))
                warning('Reference invalid – skip'); continue; end
            refRR = refMat(:,2) - refMat(1,2);     % 零基线

            try
                %% 3.1 载入模型 & 更新参数 -------------------------
                mdl = mphload(baseModel);

                % 全局 param
                mdl.param.set('x0',num2str(x0));
                mdl.param.set('y0',num2str(y0));
                mdl.param.set('a', num2str(aVal));
                mdl.param.set('force_scale',num2str(forceScale));

                % 材料 mat1 / Enu
                mat = mdl.material('mat1');
                mat.propertyGroup('Enu').set('E', sprintf('%g',E));
                mat.propertyGroup('Enu').set('nu',sprintf('%g',nu));

                % 更新 force-curve 文件 (假设函数标签仍是 int1)
                mdl.func('int1').set('filename',forceTxt{k});
                mdl.func('int1').importData;

                % 固定时间列表
                timeStep = mdl.study('std1').feature('time');
                timeStep.set('tlist',tlistStr);

                % 调宽容差（6.2 中仅 rtol）
                safeSet(mdl.sol('sol1').feature('t1'),'rtol','1e-3');

                %% 3.2 求解 --------------------------------------
                mdl.study('std1').run;
                fprintf('Solved. ');

                %% 3.3 导出 ε csv -------------------------------
                tmpCsv = fullfile(outDir,sprintf('tmp_%d.csv',k));
                if any(string(mdl.result.export.tags)=='expCsv')
                    mdl.result.export.remove('expCsv');
                end
                ex = mdl.result.export.create('expCsv','Data');
                ex.set('data','cpl1');          % 你的 cut-plane 数据集标签
                ex.set('expr',{'solid.ep1'});
                ex.set('filename',tmpCsv);
                ex.run;

                strain3D = readCsvStrain(tmpCsv); delete(tmpCsv);
                simEps   = extractStrainAt(strain3D,x0,y0,R);   % ε
                simRR    = GF*simEps - GF*simEps(1);            % ΔR/R0

                if sum(isfinite(simRR))<5
                    warning('Too few finite points – skip'); continue; end

                % 统一长度
                refRR = interp1(linspace(0,1,numel(refRR)),refRR, ...
                                linspace(0,1,numel(simRR)))';

                rmse        = sqrt(mean((simRR-refRR).^2,'omitnan'));
                rmseAll(k)  = rmse;
                fprintf('RMSE=%.3e\n',rmse);

                %% 3.4 绘图保存 -----------------------------------
                fig = figure('Visible','off');
                plot(simRR,'r','LineWidth',1.3); hold on;
                plot(refRR,'b--','LineWidth',1.3); grid on;
                xlabel('Time idx'); ylabel('ΔR/R₀ (%)');
                legend('Sim.','Ref.','Location','best');
                title(sprintf('Case %d | E=%.0e  ν=%.3f | RMSE=%.2e', ...
                      k,E,nu,rmse));
                saveas(fig,fullfile(outDir, ...
                       sprintf('cmp_%d_E%.0e_nu%.3f.png',k,E,nu)));
                close(fig);

            catch ME
                fprintf('\n❌ Case #%d failed: %s\n',k,ME.message);
            end

            % 释放模型内存
            try, ModelUtil.remove('model'); catch, end
            pause(0.05);
        end

        %% 汇总当前 (E, ν)
        avgRMSE = mean(rmseAll,'omitnan');
        fprintf('   → Avg RMSE = %.4e\n',avgRMSE);
        dlmwrite(fullfile(outDir,'summary.csv'), ...
                 [E,nu,avgRMSE],'-append','precision','%.6e');
    end
end
end
%% ===================== Helper Functions ==============================

% 若属性存在再 set；否则跳过
function safeSet(feat,prop,val)
if any(strcmp(prop,feat.properties))
    feat.set(prop,val);
end
end

% 平均给定圆邻域内的 ε
function epsVec = extractStrainAt(B,x0,y0,R)
[T,~,~] = size(B); epsVec = nan(T,1);
for t = 1:T
    x = B(t,:,1); y = B(t,:,2); v = B(t,:,3);
    idx = hypot(x-x0,y-y0)<=R;
    if any(idx), epsVec(t)=mean(v(idx),'omitnan'); end
end
end

% 解析 cut-plane CSV
function S = readCsvStrain(csvFile)
fid = fopen(csvFile,'r'); for i=1:4, fgetl(fid); end
N = sscanf(fgetl(fid),'%% Nodes,%d');
T = sscanf(fgetl(fid),'%% Expressions,%d'); for i=1:3, fgetl(fid); end
raw = zeros(N,3+T); for n = 1:N
    raw(n,:) = sscanf(fgetl(fid),'%f,').';
end
fclose(fid);
x = raw(:,1); y = raw(:,2); vs = raw(:,4:end);
S = zeros(T,N,3);
for t = 1:T
    S(t,:,1)=x.'; S(t,:,2)=y.'; S(t,:,3)=vs(:,t).';
end
end

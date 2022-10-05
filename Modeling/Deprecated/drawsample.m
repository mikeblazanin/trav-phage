function drawsample(tInd, WaveResults, SimParams, CW, MFTmodel, DrawParams, Names)
% This function draws gradient, density, speed profiles at a specified time.

% OUTPUT
% 'profiles_t_',num2str(tInd * SimParams.outDt),'.png'
% 'peak_z_t_',num2str(tInd * SimParams.outDt),'.png'
% 'c_eq_chi_df_t_',num2str(tInd*SimParams.outDt),'.png'
% 'wave_characteristics_t_',num2str(tInd*SimParams.outDt),'.png'
%% initialize
disp(['Drawing sample ', num2str(tInd), ' ...'])
z        = WaveResults.z;
dynDen   = WaveResults.dynDen;
c        = WaveResults.c;
f        = WaveResults.f;
zMax     = WaveResults.zMax;
argmaxFz = WaveResults.argmaxFz;
waveBound = DrawParams.waveBound;
%% show profiles
figure('visible', 'off')
% plot relative cell densities at selected tb
subplot(5,1,1)
plot(z, sum(dynDen)./sum(sum(dynDen)), 'LineWidth', DrawParams.lw)
hold on;
tbsShow = [0.2, 0.4, 0.5:0.01:0.6];
if SimParams.ptbType == 4
    tbsShow = tbsShow + 0.005;
end
ls = {'total'};
for tbInd = 1:numel(tbsShow)
    tbJ = tbsShow(tbInd);
    if isempty(dynDen(abs(CW.CWBias-tbJ)<0.0009,:))
        continue
    end
    plot(z, dynDen(abs(CW.CWBias-tbJ)<0.0009,:)/sum(dynDen(abs(CW.CWBias-tbJ)<0.0009,:)), 'LineWidth', DrawParams.lw);
    ls = [ls {['TB = ', num2str(tbJ)]}];
end
legend(ls)
title(['t = ', num2str(tInd * SimParams.outDt), ' s'])
xlabel('z (\mum)')
ylabel('normalized cell density')
axis tight
xlim(waveBound)

% plot absolute cell densities in a few tb bins
subplot(5,1,2)
plot(z, sum(dynDen) * SimParams.OD1)
hold on;
plot(z, sum(dynDen(abs(CW.CWBias-0.05)<0.05,:)) * SimParams.OD1, 'LineWidth', DrawParams.lw);
plot(z, sum(dynDen(abs(CW.CWBias-0.15)<0.05,:)) * SimParams.OD1, 'LineWidth', DrawParams.lw);
plot(z, sum(dynDen(abs(CW.CWBias-0.25)<0.05,:)) * SimParams.OD1, 'LineWidth', DrawParams.lw);
plot(z, sum(dynDen(abs(CW.CWBias-0.35)<0.05,:)) * SimParams.OD1, 'LineWidth', DrawParams.lw);
plot(z, sum(dynDen(abs(CW.CWBias-0.45)<0.05,:)) * SimParams.OD1, 'LineWidth', DrawParams.lw);
legend('total', 'TB = 0    ~ 0.1', 'TB = 0.1 ~ 0.2', ...
    'TB = 0.2 ~ 0.3', 'TB = 0.3 ~ 0.4', 'TB = 0.4 ~ 0.5', ...
    'location', 'northwest')
xlabel('z (\mum)')
ylabel('absolute cell density (cells/ml)')
axis tight
xlim(waveBound)

% plot aspartate profile and the perceived WaveResults.gradient
subplot(5,1,3)
plot(z, WaveResults.asp./SimParams.asp, 'LineWidth', DrawParams.lw)
hold on
plot(z, WaveResults.f./WaveResults.fMax, 'LineWidth', DrawParams.lw)
plot(z, WaveResults.gradient./WaveResults.gradientMax, 'LineWidth', DrawParams.lw)
plot(z, WaveResults.f2nd./WaveResults.f2ndMax, 'LineWidth', DrawParams.lw)
legend('s', 'f', 'df/dz', 'd^2f/dz^2', 'location', 'northwest')
xlabel('z (\mum)')
ylabel('normalized quantities')
axis tight
xlim(waveBound)

% plot potentials
subplot(5,1,4)
hold on
plot(z, c*z - MFTmodel.chi(abs(CW.CWBias - 0.05)<CW.dCWBias/2)*f, 'LineWidth', DrawParams.lw)
plot(z, c*z - MFTmodel.chi(abs(CW.CWBias - 0.1)<CW.dCWBias/2)*f,  'LineWidth', DrawParams.lw)
plot(z, c*z - MFTmodel.chi(abs(CW.CWBias - 0.15)<CW.dCWBias/2)*f, 'LineWidth', DrawParams.lw)
plot(z, c*z - MFTmodel.chi(abs(CW.CWBias - 0.2)<CW.dCWBias/2)*f,  'LineWidth', DrawParams.lw)
plot(z, c*z - MFTmodel.chi(abs(CW.CWBias - 0.25)<CW.dCWBias/2)*f, 'LineWidth', DrawParams.lw)
plot(z, c*z - MFTmodel.chi(abs(CW.CWBias - 0.3)<CW.dCWBias/2)*f,  'LineWidth', DrawParams.lw)
plot(z, c*z - MFTmodel.chi(abs(CW.CWBias - 0.35)<CW.dCWBias/2)*f, 'LineWidth', DrawParams.lw)
plot(z, c*z - MFTmodel.chi(abs(CW.CWBias - 0.4)<CW.dCWBias/2)*f,  'LineWidth', DrawParams.lw)
plot(z, c*z - MFTmodel.chi(abs(CW.CWBias - 0.45)<CW.dCWBias/2)*f, 'LineWidth', DrawParams.lw)
plot(z, c*z - MFTmodel.chi(abs(CW.CWBias - 0.5)<CW.dCWBias/2)*f,  'LineWidth', DrawParams.lw)
legend('TB = 0.05', 'TB = 0.1', 'TB = 0.15', 'TB = 0.2',...
    'TB = 0.25', 'TB = 0.3', 'TB = 0.35', 'TB = 0.4',...
    'TB = 0.45', 'TB = 0.5', 'location', 'northwest')
xlabel('z (\mum)')
ylabel('U(z)')
axis tight
xlim(waveBound)

% plot consumptions
subplot(5,1,5)
hold on
plot(z, WaveResults.oxy/SimParams.oxy, 'LineWidth', DrawParams.lw)
oxyFactor = SimParams.lAO + (1-SimParams.lAO)*1./(1+(SimParams.KAO./WaveResults.oxy).^SimParams.HAO);
plot(z, oxyFactor, 'LineWidth', DrawParams.lw)
aspOxyFactor = 1./(1+SimParams.KA./WaveResults.asp).*oxyFactor;
plot(z, aspOxyFactor, 'LineWidth', DrawParams.lw)
meanAspConsumptionFactor = sum(sum(dynDen).*aspOxyFactor) / sum(sum(dynDen));
plot(z, meanAspConsumptionFactor * ones(size(z)), '--', 'LineWidth', DrawParams.lw)
legend('oxy', 'oxy factor', 'asp oxy consumption factor', 'mean consumption factor')
xlabel('z (\mum)')
xlim(waveBound)
ylim([0 1])

set(gcf, 'PaperPosition', [0 0 1 1.1] * 3 * DrawParams.figW);
print(gcf, [Names.analysisDir, 'profiles_t_',num2str(tInd * SimParams.outDt),'.png'], '-dpng', ['-r',num2str(DrawParams.rez)], '-opengl')
close(gcf)

%% show peak position vs. chi or TB
figure('visible', 'off')
subplot(1,2,1)
hold on
plot(MFTmodel.chi, argmaxFz * ones(size(MFTmodel.chi)), 'r:', 'LineWidth', DrawParams.lw)
plot(WaveResults.chiMin * ones(size(z)), z, 'b--', 'LineWidth', DrawParams.lw)
plot(WaveResults.chiMinAtArgmaxFz * ones(size(z)), z, 'r--', 'LineWidth', DrawParams.lw)
plot(WaveResults.chiMinAtLastEps * ones(size(z)), z, 'k--', 'LineWidth', DrawParams.lw)
legend('argmax f''(z)', 'min \chi = c / max f''(z)', 'min \chi = \chi(argmax f''(z))', 'min \chi = where CDF = \epsilon')
plot(MFTmodel.chi, zMax, '.', 'LineWidth', DrawParams.lw)
title(['t = ', num2str(tInd * SimParams.outDt), ' s'])
xlabel('\chi (\mum^2/s)')
ylabel('peak z (\mum)')
axis tight
ylim(waveBound)

subplot(1,2,2)
hold on
plot(CW.CWBias, argmaxFz * ones(size(CW.CWBias)), 'r:', 'LineWidth', DrawParams.lw)
plot(WaveResults.tbMaxAtChiMin * ones(size(z)), z, 'b--', 'LineWidth', DrawParams.lw)
plot(WaveResults.tbMaxAtArgmaxFz * ones(size(z)), z, 'r--', 'LineWidth', DrawParams.lw)
plot(WaveResults.tbMaxAtLastEps * ones(size(z)), z, 'k--', 'LineWidth', DrawParams.lw)
legend('argmax f''(z)', 'max TB from min \chi = c / max f''(z)', 'max TB from min \chi = \chi(argmax f''(z))', 'max TB from min \chi = where CDF = \epsilon')
plot(CW.CWBias, zMax, '.', 'LineWidth', DrawParams.lw)
xlabel('TB')
ylabel('peak z (\mum)')
axis tight
xlim([0 1])
ylim(waveBound)

set(gcf, 'PaperPosition', [0 0 1 0.6] * DrawParams.figW * 4);
print(gcf, [Names.analysisDir, 'peak_z_t_',num2str(tInd * SimParams.outDt),'.png'], '-dpng', ['-r',num2str(DrawParams.rez)], '-opengl')
close(gcf)

%% verify c = chi f' = aA Ntot / Atot
figure('visible', 'off')
hold on;
plot(CW.CWBias, ones(size(CW.CWBias)) * c, '--', 'LineWidth', DrawParams.lw);
plot(CW.CWBias, WaveResults.gradient(WaveResults.phenotypeMaxInd)'.*MFTmodel.chi, '.', 'LineWidth', DrawParams.lw);
cInferred = meanAspConsumptionFactor * SimParams.aA * sum(sum(WaveResults.dynDen)) * SimParams.dx / SimParams.asp;
plot(CW.CWBias, ones(size(CW.CWBias)) * cInferred, '-.', 'LineWidth', DrawParams.lw);
plot(WaveResults.tbMaxAtChiMin * [1 1], [0 DrawParams.cMax], 'b--', 'LineWidth', DrawParams.lw)
plot(WaveResults.tbMaxAtArgmaxFz * [1 1], [0 DrawParams.cMax], 'r--', 'LineWidth', DrawParams.lw)
plot(WaveResults.tbMaxAtLastEps * [1 1], [0 DrawParams.cMax], 'k--', 'LineWidth', DrawParams.lw)
legend('c', '\chif''', '<\alpha_A> A_t_o_t / N_t_o_t', ...
    'max TB from min \chi = c / max f''(z)', 'max TB from min \chi = \chi(argmax f''(z))', 'max TB from min \chi = where CDF = \epsilon')
title(['t = ', num2str(tInd*SimParams.outDt), ' s'])
xlabel('TB')
ylabel('speed (\mum/s)')
xlim([0 1])
ylim([0 DrawParams.cMax])

set(gcf, 'PaperPosition', [0 0 1 1] * DrawParams.figW);
print(gcf, [Names.analysisDir, 'c_eq_chi_df_t_',num2str(tInd*SimParams.outDt),'.png'], '-dpng', ['-r',num2str(DrawParams.rez)], '-opengl')
close(gcf)

%% show density, gradient, and speed in moving coorinate
nPheno = sum(dynDen, 2);
PWaveT = nPheno / sum(nPheno) / CW.dCWBias;
% PWaveTExcluded = PWaveT;
% PWaveTExcluded(WaveResults.tbMaxAtArgmaxFz/CW.dCWBias:end) = 0;
% PWaveTExcluded = PWaveTExcluded/sum(PWaveTExcluded)/CW.dCWBias;
populations = {};
populations{1} = dynDen(abs(CW.CWBias-0.1)<0.1,:);
populations{2} = dynDen(abs(CW.CWBias-0.25)<0.05,:);
populations{3} = dynDen(abs(CW.CWBias-0.35)<0.05,:);
populations{4} = dynDen(abs(CW.CWBias-0.45)<0.05,:);
zPopInd = {};
zMaxPopInd = {};
for i = 1:4
    [~, zPopInd{i}] = max(sum(populations{i}));
    [~, zMaxPopInd{i}] = min(abs(zMax - z(zPopInd{i})));
end
colors = {'r', 'y', 'g', 'c'};

figure('visible', 'off')
subplot(5,1,1)
hold on
plot(CW.CWBias, PWaveT, 'LineWidth', DrawParams.lw)
plot(WaveResults.tbMaxAtChiMin*[1 1], [0 10], 'b--', 'LineWidth', DrawParams.lw)
plot(WaveResults.tbMaxAtArgmaxFz*[1 1], [0 10], 'r--', 'LineWidth', DrawParams.lw)
plot(WaveResults.tbMaxAtLastEps*[1 1], [0 10], 'k--', 'LineWidth', DrawParams.lw)
legend('max TB from min \chi = c / max f''(z)', 'max TB from min \chi = \chi(argmax f''(z))', 'max TB from min \chi = where CDF = \epsilon')
xlim([0 1])
ylim([0 10])
xlabel('Tumble bias')
ylabel('P(TB)')

subplot(5,1,2)
hold on
exp1 = MFTmodel.chi .* PWaveT * CW.dCWBias ./ abs([diff(MFTmodel.chi); 0]);
exp2 = ones(size(CW.CWBias)) * SimParams.KiA/SimParams.asp;
exp3 = exp2 + cumsum(PWaveT,'reverse') * CW.dCWBias;
plot(CW.CWBias, exp1, 'LineWidth', DrawParams.lw)
plot(CW.CWBias, exp2, '--', 'LineWidth', DrawParams.lw)
plot(CW.CWBias, exp3, 'LineWidth', DrawParams.lw)
plot(CW.CWBias, exp3./exp1, 'LineWidth', DrawParams.lw)
plot(WaveResults.tbMaxAtChiMin*[1 1], [0 1.5], 'b--', 'LineWidth', DrawParams.lw)
plot(WaveResults.tbMaxAtArgmaxFz*[1 1], [0 1.5], 'r--', 'LineWidth', DrawParams.lw)
plot(WaveResults.tbMaxAtLastEps*[1 1], [0 1.5], 'k--', 'LineWidth', DrawParams.lw)
legend('(1)=\chiP(\chi)', '(2)=K_i/A_\infty', '(3)=(2)+\int_0^\chiP(\chi'')d\chi'' ', '(3)/(1)', ...
    'max TB from min \chi = c / max f''(z)', 'max TB from argmax f''(z)', 'max TB from min \chi = where CDF = \epsilon')
xlim([0 1])
ylim([0 1.5])
xlabel('Tumble bias')

subplot(5,1,3)
plot(z, sum(dynDen) * SimParams.OD1 / 1e9, 'b')
hold on;
for i = 1:4
    plot(z, sum(populations{i}) * SimParams.OD1 / 1e9, '--', 'color', colors{i}, 'LineWidth', DrawParams.lw);
end
plot(z, log((1+WaveResults.asp/SimParams.KiA)./(1+WaveResults.asp/SimParams.KaA)), 'k', 'LineWidth', DrawParams.lw)
plot(argmaxFz*[1 1], [0 8], 'r:', 'LineWidth', DrawParams.lw)
xlabel('z (\mum)')
ylabel('\rho (10^9 ml^-^1), f(z)')
xlim(waveBound)
ylim([0 15])

subplot(5,1,4)
hold on
plot(zMax, MFTmodel.chi/3000, 'b', 'LineWidth', DrawParams.lw)
plot(z, WaveResults.gradient*1000, 'k', 'LineWidth', DrawParams.lw)
for i = 1:4
    scatter(zMax(zMaxPopInd{i}), MFTmodel.chi(zMaxPopInd{i})/3000, 50, colors{i}, 'filled')
    scatter(z(zPopInd{i}), WaveResults.gradient(zPopInd{i})*1000, 50, colors{i}, 'filled')
end
plot(argmaxFz*[1 1], [0 3], 'r:', 'LineWidth', DrawParams.lw)
xlabel('z (\mum)')
ylabel('\chi (3x10^3 \mum^2/s), df/dz (mm^-^1)')
xlim(waveBound)

subplot(5,1,5)
hold on
plot(zMax, WaveResults.gradient(WaveResults.phenotypeMaxInd)'.*MFTmodel.chi, 'k', 'LineWidth', DrawParams.lw);
plot(z, ones(size(z)).*c, 'k--', 'LineWidth', DrawParams.lw);
for i = 1:4
    scatter(zMax(zMaxPopInd{i}), WaveResults.gradient(WaveResults.phenotypeMaxInd(zMaxPopInd{i})).*MFTmodel.chi(zMaxPopInd{i}), 50, colors{i}, 'filled')
end
plot(argmaxFz*[1 1], [0 4], 'r:', 'LineWidth', DrawParams.lw)
ylabel('c (\mum/s)')
xlabel('z (\mum)')
xlim(waveBound)
ylim([0 10])

set(gcf, 'PaperPosition', [0 0 1 3] * DrawParams.figW);
print(gcf, [Names.analysisDir, 'wave_characteristics_t_',num2str(tInd*SimParams.outDt),'.png'], '-dpng', ['-r',num2str(DrawParams.rez)], '-opengl')
close(gcf)
end
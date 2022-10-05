function drawphenoprofilemovieframe(tInd, WaveResults, ...
    SimResults, SimParams, CW, MFTmodel, DrawParams, Names)
% This function saves image frames of the profiles of selected phenotypes
% at a certain time.

disp(['Saving phenotype profiles movie frame ', num2str(tInd), ' ...'])
showInds = 200:50:700;
if SimParams.ptbType == 4
    showInds = 10:10:140;
end
figure('visible', 'off')

% plot the density profiles
hold on
plot(SimParams.x*1000, squeeze( SimResults.dynProfile.cellDensity(showInds,:,tInd) )'*SimParams.OD1, 'LineWidth', DrawParams.lw)
% plot(SimParams.x*1000, squeeze( sum(SimResults.dynProfile.cell_density(:,:,t)) )*SimParams.OD1, 'k--', 'LineWidth', DrawParams.lw)
ls = [];
for tbInd = 1:numel(showInds)
    tbJ = CW.CWBias(showInds(tbInd));
    ls = [ls {['TB = ', num2str(tbJ)]}];
end
% ls = [ls {'Total'}];
legend(ls)
xlabel('x (mm)')
ylabel('cell density (cells/ml)')
title(['t = ', num2str(tInd*SimParams.outDt/60, '%.1f'), ' min'])
axis tight
ylim([0 0.05*SimParams.OD1])

% plot an inset showing c = chi * f'
axes('Position', [0.55 0.65 0.25 0.25])
box on
plot(CW.CWBias, WaveResults.gradient(WaveResults.phenotypeMaxInd)'.*MFTmodel.chi, '.', 'LineWidth', DrawParams.lw);
hold on;
plot(CW.CWBias, ones(size(CW.CWBias)).*WaveResults.c, '--', 'LineWidth', DrawParams.lw);
plot(WaveResults.tbMaxAtChiMin * [1 1], [0 DrawParams.cMax], 'k--', 'LineWidth', DrawParams.lw)
plot(WaveResults.tbMaxAtArgmaxFz * [1 1], [0 DrawParams.cMax], 'g-.', 'LineWidth', DrawParams.lw)
legend('\chif''', 'c')
xlabel('TB')
ylabel('speed (\mum/s)')
ylim([0 DrawParams.cMax])

% save plot
set(gcf, 'PaperPosition', [0 0 1 0.6] * DrawParams.figW * 4);
if ~exist(Names.phenoProfileMovieDir, 'dir')
    mkdir(Names.phenoProfileMovieDir)
end
print(gcf, [Names.phenoProfileMovieDir, 'pheno_profile_',num2str(tInd*SimParams.outDt),'.jpg'], '-djpeg', ['-r',num2str(DrawParams.rez)], '-opengl')
close(gcf)
end
function drawleakagedynamics(NumberDensities, SimParams, CW, DrawParams, Names)
% This function draws simulated and predicted number densitites over time.

% OUTPUT
% leakagedynamics.png

disp('Drawing leakage dynamics ...')
saveDir = [Names.analysisDir, 'leakage\'];
if ~exist(saveDir, 'dir')
    mkdir(saveDir)
end

figure('visible', 'off')
hold on
dfc = get(gca,'colororder'); % default colors
nColors = size(dfc, 1);
i = 1;
for t_ind = DrawParams.tInds
    plot(CW.CWBias, NumberDensities.simulatedND(:,t_ind), 'Color', dfc(i,:), 'LineWidth', DrawParams.lw/3)
    i = mod(i, nColors) + 1;
end
i = 1;
for t_ind = DrawParams.tInds
    plot(CW.CWBias, NumberDensities.predictedND(:,t_ind), 'Color', dfc(i,:), 'LineStyle', '--', 'LineWidth', DrawParams.lw)
    i = mod(i, nColors) + 1;
end
xlabel('TB')
ylabel('pdf')
axis tight
xlim([0 1])
ls = {};
for t_ind = 1:DrawParams.nTInds
    ls = [ls {['t = ', num2str(DrawParams.tInds(t_ind) * SimParams.outDt), ' s']}];
end
legend(ls)

set(gcf, 'PaperPosition', [0 0 1 0.7] * 3 * DrawParams.figW);
figName = [saveDir, 'leakagedynamics.png'];
disp(['Saving figure ', figName, '...'])
print(gcf, figName, '-dpng', ['-r',num2str(DrawParams.rez)], '-opengl')
close(gcf)
end
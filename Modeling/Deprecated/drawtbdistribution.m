function drawtbdistribution(AnalysisResults, SimParams, CW, DrawParams, Names)
% This function draws tumble bias distribution and fold enrichment over
% time.

% OUTPUT
% tb_dist.png

disp('Drawing tumble bias distributions ...')
figure('visible', 'off')
subplot(1,2,1)
plot(CW.CWBias, AnalysisResults.PWaveAll(:,DrawParams.tInds), 'LineWidth', DrawParams.lw)
hold on
plot(CW.CWBias, CW.P, '--k', 'LineWidth', DrawParams.lw)
xlabel('TB')
ylabel('pdf')
axis tight
xlim([0 1])
ls = {};
for tInd = 1:DrawParams.nTInds
    ls = [ls {['t = ', num2str(DrawParams.tInds(tInd) * SimParams.outDt), ' s']}];
end
ls = [ls {'t = 0'}];
legend(ls)

subplot(1,2,2)
plot(CW.CWBias, AnalysisResults.PWaveAll(:,DrawParams.tInds)./repmat(CW.P,1,DrawParams.nTInds), 'LineWidth', DrawParams.lw)
xlabel('TB')
ylabel('Fold Enrichment')
axis tight
xlim([0 1])
limsy=get(gca,'YLim');
set(gca,'Ylim',[0 limsy(2)]);


set(gcf, 'PaperPosition', [0 0 1 0.3] * 3 * DrawParams.figW);
figName = [Names.analysisDir, 'tb_dist.png'];
disp(['Saving figure ', figName, '...'])
print(gcf, figName, '-dpng', ['-r',num2str(DrawParams.rez)], '-opengl')
close(gcf)
end
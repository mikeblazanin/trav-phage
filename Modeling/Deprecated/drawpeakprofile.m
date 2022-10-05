function drawpeakprofile(AnalysisResults, SimParams, DrawParams, Names,  ...
    plotMultiple, nCases, AnalysisResultsMultiple)
% This function draws temporal evolution of wave peaks.

% OUTPUT
% analysis.png

disp('Drawing peak profiles ...')
if plotMultiple
%     lg = {'Asp 50', 'Asp 100', 'Asp 200'};
%     lg = {'Hill Conc Bell TB', 'Hill Conc Unif TB', 'Flat Conc Bell TB', 'Flat Conc Unif TB'};
%     lg = {'dx = 10 dt = 0.02', 'dx = 20 dt = 0.1', 'dx = 25 dt = 0.125', 'dx = 50 dt = 0.5'};
%     lg = {'TB mean 0.5 std 0.1', 'TB mean 0.3 std 0.26', 'TB mean 0.3 std 0.01', 'TB mean 0.3 std 0.01 2x cells'};
%     lg = {'TB mean 0.5 std 0.1', 'TB mean 0.3 std 0.1', 'TB mean 0.3 std 0.1 Ki=0.5', 'TB mean 0.3 std 0.01', 'TB mean 0.3 std 0.01 2x cells', 'TB mean 0.3 std 0.01 Ki=0.5'};
    lg = {'P(TB) wide', 'P(TB) narrow'};
end

figure('visible', 'off')
subplot(2,2,1)
if plotMultiple
    hold on
    for ncase = 1:nCases
        plot(SimParams.time/60, AnalysisResultsMultiple.peakSpdMultiple{ncase}, 'LineWidth', DrawParams.lw)
    end
else
    plot(SimParams.time/60, AnalysisResults.peakSpd, 'LineWidth', DrawParams.lw)
end
title('peak speed')
xlabel('time (min)')
ylabel('speed (\mum/s)')
axis tight
ylim([0 15])

subplot(2,2,2)
if plotMultiple
    hold on
    for ncase = 1:nCases
        plot(SimParams.time/60, AnalysisResultsMultiple.peakDenMultiple{ncase}, 'LineWidth', DrawParams.lw)
    end
    legend(lg)
else
    plot(SimParams.time/60, AnalysisResults.peakDen, 'LineWidth', DrawParams.lw)
end
title('peak density')
xlabel('time (min)')
ylabel('density (cells/ml)')
axis tight
limsy=get(gca,'YLim');
set(gca,'Ylim',[0 limsy(2)]);

subplot(2,2,3)
if plotMultiple
    hold on
    for ncase = 1:nCases
        plot(SimParams.time/60, AnalysisResultsMultiple.tbExcludeMultiple{ncase}, 'LineWidth', DrawParams.lw)
    end
else
    plot(SimParams.time/60, AnalysisResults.tbExclude, 'LineWidth', DrawParams.lw)
end
title('boundary tumble bias')
xlabel('SimParams.time (min)')
ylabel('tumble bias')
axis tight
ylim([0 1])

subplot(2,2,4)
if plotMultiple
    hold on
    for ncase = 1:nCases
        plot(SimParams.time/60, AnalysisResultsMultiple.waveCellsMultiple{ncase}, 'LineWidth', DrawParams.lw)
    end
else
    plot(SimParams.time/60, AnalysisResults.waveCellsAll, 'LineWidth', DrawParams.lw)
end
title('total cells within wave')
xlabel('time (min)')
ylabel('cells')
axis tight
limsy=get(gca,'YLim');
set(gca,'Ylim',[0 limsy(2)]);

set(gcf, 'PaperPosition', [0 0 1 0.6] * 2 * DrawParams.figW);
figName = [Names.analysisDir, 'analysis.png'];
disp(['Saving figure ', figName, '...'])
print(gcf, figName, '-dpng', ['-r',num2str(DrawParams.rez)], '-opengl')
close(gcf)
end
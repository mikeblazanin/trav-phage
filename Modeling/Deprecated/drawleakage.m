function drawleakage(AnalysisParams, AnalysisResults, SimParams, CW, MFTmodel, DrawParams, Names, plotMultiple)
% This function draws cell leakage using various methods.

% OUTPUT
% leakage.png
% relative_leakage.png
% 'leakage_t=', num2str(tInd * SimParams.outDt), '.png'

% initialization
disp('Drawing cell leakage ...')
thinning = 1;
ncwbThinned = floor((CW.ncwb-1)/thinning)+1;
x3 = repmat(SimParams.time, ncwbThinned,1)';
y3 = repmat(CW.CWBias(1:thinning:end), 1, SimParams.nT)';
waveCellsPhenoThinned = repmat(AnalysisResults.waveCellsGMax, ncwbThinned,1) .* AnalysisResults.PWaveGMax(1:thinning:end,:) * CW.dCWBias;
z31 = waveCellsPhenoThinned';
z32 = -[zeros(ncwbThinned,AnalysisParams.calculatePeakSpacing) waveCellsPhenoThinned(:,AnalysisParams.calculatePeakSpacing+1:end)-waveCellsPhenoThinned(:,1:end-AnalysisParams.calculatePeakSpacing)]'/AnalysisParams.calculatePeakSpacing/SimParams.outDt;
z33 = AnalysisResults.fluxBackCells(1:thinning:end,:)';
x3 = x3(10:end,:);
y3 = y3(10:end,:);
z31 = z31(10:end,:);
z32 = z32(10:end,:);
z33 = z33(10:end,:);

waveCellsPhenoGMax = repmat(AnalysisResults.waveCellsGMax, CW.ncwb,1) .* AnalysisResults.PWaveGMax * CW.dCWBias;
fluxBackSim = -[zeros(CW.ncwb,AnalysisParams.calculatePeakSpacing) waveCellsPhenoGMax(:,AnalysisParams.calculatePeakSpacing+1:end)-waveCellsPhenoGMax(:,1:end-AnalysisParams.calculatePeakSpacing)]/AnalysisParams.calculatePeakSpacing/SimParams.outDt;

saveDir = [Names.analysisDir, 'leakage\'];
if ~exist(saveDir, 'dir')
    mkdir(saveDir)
end

figure('visible', 'off')
subplot(3,3,1)
plot3(x3, y3, z31, 'LineWidth', DrawParams.lw)
xlim([0 SimParams.nT*SimParams.outDt])
ylim([0 0.6])
xlabel('time (s)')
ylabel('TB')
zlabel('cell number (sim)')
subplot(3,3,2)
surf(x3, y3, z31, 'LineWidth', DrawParams.lw, 'EdgeColor','none')
%         colorbar
xlim([0 SimParams.nT*SimParams.outDt])
ylim([0 0.6])
xlabel('time (s)')
ylabel('TB')
zlabel('cell number (sim)')
subplot(3,3,3)
plot(SimParams.time(10:end), sum(z31,2), 'LineWidth', DrawParams.lw)
xlim([0 SimParams.nT*SimParams.outDt])
ylim([0 SimParams.NCells])
xlabel('time (s)')
ylabel('cell number (sim)')

subplot(3,3,4)
plot3(x3, y3, z32, 'LineWidth', DrawParams.lw)
xlim([0 SimParams.nT*SimParams.outDt])
ylim([0 0.6])
zlim([-0.1 2])
xlabel('time (s)')
ylabel('TB')
zlabel('number flux (/s) (sim)')

subplot(3,3,5)
surf(x3, y3, z32, 'LineWidth', DrawParams.lw, 'EdgeColor','none')
%         colorbar
caxis([-0.1 2])
xlim([0 SimParams.nT*SimParams.outDt])
ylim([0 0.6])
zlim([-0.1 2])
xlabel('time (s)')
ylabel('TB')
zlabel('number flux (/s) (sim)')

subplot(3,3,6)
plot(SimParams.time(10:end), sum(z32,2), 'LineWidth', DrawParams.lw)
xlim([0 SimParams.nT*SimParams.outDt])
ylim([-0.1 20])
xlabel('time (s)')
ylabel('number flux (/s) (sim)')

subplot(3,3,7)
plot3(x3, y3, z33, 'LineWidth', DrawParams.lw)
xlim([0 SimParams.nT*SimParams.outDt])
ylim([0 0.6])
zlim([-0.1 2])
xlabel('time (s)')
ylabel('TB')
zlabel('number flux (/s) (theor)')

subplot(3,3,8)
surf(x3, y3, z33, 'LineWidth', DrawParams.lw, 'EdgeColor','none')
%         colorbar
caxis([-0.1 2])
xlim([0 SimParams.nT*SimParams.outDt])
ylim([0 0.6])
zlim([-0.1 2])
xlabel('time (s)')
ylabel('TB')
zlabel('number flux (/s) (theor)')

subplot(3,3,9)
plot(SimParams.time(10:end), sum(z33,2), 'LineWidth', DrawParams.lw)
xlim([0 SimParams.nT*SimParams.outDt])
ylim([-0.1 20])
xlabel('time (s)')
ylabel('number flux (/s) (theor)')

set(gcf, 'PaperPosition', [0 0 1 0.6] * 4 * DrawParams.figW);
figName = [saveDir, 'leakage.png'];
disp(['Saving figure ', figName, '...'])
print(gcf, figName, '-dpng', ['-r',num2str(DrawParams.rez)], '-opengl')
close(gcf)

% plot relative leakage
totalBackFlux = sum(fluxBackSim);
bdAtLastEpsChi = (AnalysisResults.chiMinAtLastEps-AnalysisResults.muMinAtLastEps).^2 ./ AnalysisResults.chiMinAtLastEps;
bdagmEqChi = (AnalysisResults.chiMinAtArgmaxFz-AnalysisResults.muMinAtArgmaxFz).^2 ./ AnalysisResults.chiMinAtArgmaxFz;
bdcgmEqChi = (AnalysisResults.chiMin-AnalysisResults.muMinAtChiMin).^2 ./ AnalysisResults.chiMin;
lkEqChi = ((MFTmodel.chi-MFTmodel.mu).^2./MFTmodel.chi)' * fluxBackSim ./ totalBackFlux;
meanEqChi = ((MFTmodel.chi-MFTmodel.mu).^2./MFTmodel.chi)' * AnalysisResults.PWaveGMax * CW.dCWBias;
bdAtLastEpsTimeScale = bdAtLastEpsChi ./ AnalysisResults.peakSpd.^2;
bdagmTimeScale = bdagmEqChi ./ AnalysisResults.peakSpd.^2;
bdcgmTimeScale = bdcgmEqChi ./ AnalysisResults.peakSpd.^2;
lkTimeScale = lkEqChi ./ AnalysisResults.peakSpd.^2;
meanTimeScale = meanEqChi ./ AnalysisResults.peakSpd.^2;
meanRateScale = AnalysisResults.waveCellsGMax ./ meanTimeScale;

if plotMultiple
    figure(100)
    lg = {'TB mean 0.5 std 0.1', 'TB mean 0.3 std 0.1', 'TB mean 0.3 std 0.1 Ki=0.5', 'TB mean 0.3 std 0.01', 'TB mean 0.3 std 0.01 2x cells', 'TB mean 0.3 std 0.01 Ki=0.5'};
else
    figure('visible', 'off')
end
subplot(4,1,1)
if plotMultiple
    hold on
end
plot(SimParams.time, AnalysisResults.waveCellsGMax, 'LineWidth', DrawParams.lw)
xlim([0 SimParams.nT*SimParams.outDt])
ylim([0 SimParams.NCells])
xlabel('time (s)')
ylabel('cell number')
if plotMultiple
    legend(lg)
end

subplot(4,1,2)
if plotMultiple
    hold on
end
plot(SimParams.time, totalBackFlux, 'LineWidth', DrawParams.lw)
xlim([0 SimParams.nT*SimParams.outDt])
ylim([0 20])
xlabel('time (s)')
ylabel('number flux (/s)')

subplot(4,1,3)
if plotMultiple
    hold on
end
hold on
plot(SimParams.time, totalBackFlux./AnalysisResults.waveCellsGMax, 'LineWidth', DrawParams.lw)
plot(SimParams.time, SimParams.KiA/SimParams.asp/2./bdAtLastEpsTimeScale, 'k', 'LineWidth', DrawParams.lw)
plot(SimParams.time, SimParams.KiA/SimParams.asp/2./bdagmTimeScale, 'y', 'LineWidth', DrawParams.lw)
plot(SimParams.time, SimParams.KiA/SimParams.asp/2./bdcgmTimeScale, 'c', 'LineWidth', DrawParams.lw)
plot(SimParams.time, SimParams.KiA/SimParams.asp/2./lkTimeScale, 'r', 'LineWidth', DrawParams.lw)
plot(SimParams.time, SimParams.KiA/SimParams.asp/2./meanTimeScale, 'g', 'LineWidth', DrawParams.lw)
legend('|dlog(N)/dt|', ...
    '\epsilon/2/T(\chi = where CDF = \epsilon)', ...
    '\epsilon/2/T(\chi peaked at gradient max)', ...
    '\epsilon/2/T(\chi = c / gradient max)', ...
    '\epsilon/2/T(\chi = weighted mean of leaking cells)', ...
    '\epsilon/2/T(\chi = weighted mean of all cells)')
xlim([0 SimParams.nT*SimParams.outDt])
ylim([0 1] * 1e-3)
xlabel('time (s)')
ylabel('relative flux (/s)')

subplot(4,1,4)
if plotMultiple
    hold on
end
plot(SimParams.time, totalBackFlux./meanRateScale, 'LineWidth', DrawParams.lw)
xlim([0 SimParams.nT*SimParams.outDt])
ylim([0 0.05])
xlabel('time (s)')
ylabel('non-dimensional flux')

if ~plotMultiple
    set(gcf, 'PaperPosition', [0 0 1 0.8] * 4 * DrawParams.figW);
    figName = [saveDir, 'relative_leakage.png'];
    disp(['Saving figure ', figName, '...'])
    print(gcf, figName, '-dpng', ['-r',num2str(DrawParams.rez)], '-opengl')
    close(gcf)
end

% plot TB specific leakage over snapshots
for tInd = DrawParams.tInds
    figure('visible', 'off')
    hold on
    plot(CW.CWBias, fluxBackSim(:,tInd), 'k', 'LineWidth', DrawParams.lw)
    plot(CW.CWBias, -AnalysisResults.rhoLinMaxGrad(:,tInd) .* AnalysisResults.relSpeedMaxGrad(:,tInd), 'k--', 'LineWidth', DrawParams.lw)
    plot(CW.CWBias, AnalysisResults.fluxBackCells(:,tInd), 'm', 'LineWidth', DrawParams.lw)
    plot(CW.CWBias, AnalysisResults.imFluxBackCells(:,tInd), 'm--', 'LineWidth', DrawParams.lw)
    plot(CW.CWBias, AnalysisResults.acFluxBackCells(:,tInd), 'r:', 'LineWidth', DrawParams.lw)
    plot(CW.CWBias, AnalysisResults.gsFluxBackCells(:,tInd,1), 'y', 'LineWidth', DrawParams.lw)
    plot(CW.CWBias, AnalysisResults.gsFluxBackCells(:,tInd,2), 'c', 'LineWidth', DrawParams.lw)
%     plot(CW.CWBias, AnalysisResults.gsFluxBackCells(:,tInd,3), 'g', 'LineWidth', DrawParams.lw)
%     plot(CW.CWBias, AnalysisResults.gsFluxBackCells(:,tInd,4), 'r', 'LineWidth', DrawParams.lw)
%     plot(CW.CWBias, AnalysisResults.gsFluxBackCells(:,tInd,5), 'g', 'LineWidth', DrawParams.lw)
    plot(CW.CWBias, AnalysisResults.gsFluxBackCells(:,tInd,6), 'g', 'LineWidth', DrawParams.lw)
    plot(CW.CWBias, AnalysisResults.gsFluxBackCells(:,tInd,7), 'b', 'LineWidth', DrawParams.lw)
%     plot(CW.CWBias, AnalysisResults.cpFluxBackCells(:,tInd), 'b')
%     plot(CW.CWBias, AnalysisResults.cqFluxBackCells(:,tInd), 'b--')
%     plot(CW.CWBias, AnalysisResults.fpFluxBackCells(:,tInd), 'g')
%     plot(CW.CWBias, AnalysisResults.fqFluxBackCells(:,tInd), 'g--')
%     if tInd-1 > 0
%         plot(CW.CWBias, AnalysisResults.fluxBackCells(:,tInd-1), 'r--')
%         plot(CW.CWBias, AnalysisResults.imFluxBackCells(:,tInd-1), 'r:')
%         plot(CW.CWBias, AnalysisResults.acFluxBackCells(:,tInd-1), 'r-.')
%     end
%     if tInd+1 <= SimParams.nT
%         plot(CW.CWBias, AnalysisResults.fluxBackCells(:,tInd+1), 'b--')
%         plot(CW.CWBias, AnalysisResults.imFluxBackCells(:,tInd+1), 'b:')
%         plot(CW.CWBias, AnalysisResults.acFluxBackCells(:,tInd+1), 'b-.')
%     end
%     legend('simulation','theory (Kramers)','theory (Kramers higher order)', ...
%         'theory (non-Gaussian potential)', ...
%         'c\rho(group boundary)', 'c\rho(individual boundary)', ...
%         '\chif''\rho(group boundary)', '\chif''\rho(individual boundary)')
    legend([num2str(nansum(fluxBackSim(:,tInd)),'%.2e'),' simulation (integration)'], ...
        [num2str(nansum(-AnalysisResults.rhoLinMaxGrad(:,tInd) .* AnalysisResults.relSpeedMaxGrad(:,tInd)),'%.2e'),' simulation (flux)'], ...
        [num2str(nansum(AnalysisResults.fluxBackCells(:,tInd)),'%.2e'),' theory (Kramers)'], ...
        [num2str(nansum(AnalysisResults.imFluxBackCells(:,tInd)),'%.2e'),' theory (Kramers higher order)'], ...
        [num2str(nansum(AnalysisResults.acFluxBackCells(:,tInd)),'%.2e'),' theory (non-Gaussian potential)'], ...
        [num2str(nansum(AnalysisResults.gsFluxBackCells(:,tInd,1)),'%.2e'),' theory (\chi_b peaked at gradient max)'], ...
        [num2str(nansum(AnalysisResults.gsFluxBackCells(:,tInd,2)),'%.2e'),' theory (\chi_b = c / gradient max)'], ...
        [num2str(nansum(AnalysisResults.gsFluxBackCells(:,tInd,6)),'%.2e'),' theory (\chi_b at CDF = \epsilon)'], ...
        [num2str(nansum(AnalysisResults.gsFluxBackCells(:,tInd,7)),'%.2e'),' theory (\chi_b = where c = speed limit)'])
%         [num2str(nansum(AnalysisResults.gsFluxBackCells(:,tInd,6)),'%.2e'),' theory (\chi_b at |d^2f/dz^2| << (df/dz)^2)'], ...
%         'theory (\chi_b at A >> K_i)', ...
%         'theory (\chi_b at A >> sqrt(K_iK_A))', ...
    xlabel('TB')
    ylabel('number flux (/s)')
    xlim([0 0.7])
    ylim([-0.1 2])
    title(['t = ', num2str(tInd * SimParams.outDt), ' s'])
    
    set(gcf, 'PaperPosition', [0 0 1 0.6] * 2 * DrawParams.figW);
    figName = [saveDir, 'leakage_t=', num2str(tInd * SimParams.outDt), '.png'];
    disp(['Saving figure ', figName, '...'])
    print(gcf, figName, '-dpng', ['-r',num2str(DrawParams.rez)], '-opengl')
    close(gcf)
end
end
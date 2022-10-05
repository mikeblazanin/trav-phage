function drawgradmaxdrift(AnalysisResults, SimParams, CW, DrawParams, Names)
% This function draws drift of phenotypes at maximal gradient.

% OUTPUT
% 'grad_max_drift_t=', num2str(t * SimParams.outDt), '.png'
% 'grad_max_drift_TB=', num2str(CW.CWBias(cwInd),'%.3f'), '.png'

disp('Drawing drifts at gradient maximum ...')
saveDir = [Names.analysisDir, 'grad_max_drift\'];
if ~exist(saveDir, 'dir')
    mkdir(saveDir)
end

for t = DrawParams.tInds
    figure('visible', 'off')
    subplot(4,1,1)
    hold on
    plot(CW.CWBias, ones(CW.ncwb,1) * AnalysisResults.dxdtMaxGrad(t), 'r', 'LineWidth', DrawParams.lw)
    plot(CW.CWBias, AnalysisResults.cfDriftMaxGrad(:,t), 'b', 'LineWidth', DrawParams.lw)
    plot(CW.CWBias, AnalysisResults.cfDriftPhePeak(:,t), 'b--', 'LineWidth', DrawParams.lw)
    plot(CW.CWBias, AnalysisResults.mrDriftMaxGrad(:,t), 'm', 'LineWidth', DrawParams.lw)
    plot(CW.CWBias, AnalysisResults.relSpeedMaxGrad(:,t), 'k', 'LineWidth', DrawParams.lw)
    plot(CW.CWBias, zeros(CW.ncwb,1), 'k--', 'LineWidth', DrawParams.lw)
    legend('dxL/dt', '\chi_idf/dz', '\chi_idf/dz(z_i)', '-\mu_idln(\rho_i)/dz', 'v_i')
    xlabel('TB')
    ylabel('drift speed (\mum/s)')
    xlim([0 0.7])
    ylim([-5 20])
    title(['t = ', num2str(t * SimParams.outDt), ' s'])
    
    subplot(4,1,2)
    hold on
    plot(CW.CWBias, AnalysisResults.peakPhenoPos(:,t) - AnalysisResults.xMaxGrad(t), 'r', 'LineWidth', DrawParams.lw)
    plot(CW.CWBias, zeros(size(CW.CWBias)), 'r--', 'LineWidth', DrawParams.lw)
    xlabel('TB')
    ylabel('peak position z_i (\mum)')
    xlim([0 0.7])
    ylim([-2 3] * 1e3)
    
    subplot(4,1,3)
    hold on
    plot(CW.CWBias, AnalysisResults.rhoLinMaxGrad(:,t) .* ones(CW.ncwb,1) * AnalysisResults.dxdtMaxGrad(t), 'r', 'LineWidth', DrawParams.lw)
    plot(CW.CWBias, AnalysisResults.rhoLinMaxGrad(:,t) .* AnalysisResults.cfDriftMaxGrad(:,t), 'b', 'LineWidth', DrawParams.lw)
    plot(CW.CWBias, AnalysisResults.rhoLinPhePeak(:,t) .* AnalysisResults.cfDriftPhePeak(:,t), 'b--', 'LineWidth', DrawParams.lw)
    plot(CW.CWBias, AnalysisResults.rhoLinMaxGrad(:,t) .* AnalysisResults.mrDriftMaxGrad(:,t), 'm', 'LineWidth', DrawParams.lw)
    plot(CW.CWBias, AnalysisResults.rhoLinMaxGrad(:,t) .* AnalysisResults.relSpeedMaxGrad(:,t), 'k', 'LineWidth', DrawParams.lw)
    legend('\rho_idxL/dt', '\rho_i\chi_idf/dz', '\rho_i\chi_idf/dz(z_i)', '-\mu_id\rho_i/dz', '\rho_iv_i')
    xlabel('TB')
    ylabel('number flux (/s)')
    xlim([0 0.7])
    ylim([-0.1 1] * 1e2)
    
    subplot(4,1,4)
    hold on
    plot(CW.CWBias, AnalysisResults.cfDivMaxGrad(:,t), 'b', 'LineWidth', DrawParams.lw)
    plot(CW.CWBias, AnalysisResults.mrDivMaxGrad(:,t), 'm', 'LineWidth', DrawParams.lw)
    plot(CW.CWBias, AnalysisResults.cfDivMaxGrad(:,t) + AnalysisResults.mrDivMaxGrad(:,t), 'k', 'LineWidth', DrawParams.lw)
    legend('\chi_id^2f/dz^2', '-\mu_id^2\rho_i/dz^2', 'dv_i/dz', 'location', 'SouthEast')
    xlabel('TB')
    ylabel('velocity divergence (/s)')
    xlim([0 0.7])
    ylim([-8 3] * 1e-3)
    
    set(gcf, 'PaperPosition', [0 0 2 4] * DrawParams.figW);
    figName = [saveDir, 'grad_max_drift_t=', num2str(t * SimParams.outDt), '.png'];
    disp(['Saving figure ', figName, '...'])
    print(gcf, figName, '-dpng', ['-r',num2str(DrawParams.rez)], '-opengl')
    close(gcf)
end

for cwInd = 20:1:40
    figure('visible', 'off')
    subplot(4,1,1)
    hold on
    plot(SimParams.time, AnalysisResults.dxdtMaxGrad, 'r', 'LineWidth', DrawParams.lw)
    plot(SimParams.time, AnalysisResults.cfDriftMaxGrad(cwInd,:), 'b', 'LineWidth', DrawParams.lw)
    plot(SimParams.time, AnalysisResults.cfDriftPhePeak(cwInd,:), 'b--', 'LineWidth', DrawParams.lw)
    plot(SimParams.time, AnalysisResults.mrDriftMaxGrad(cwInd,:), 'm', 'LineWidth', DrawParams.lw)
    plot(SimParams.time, AnalysisResults.relSpeedMaxGrad(cwInd,:), 'k', 'LineWidth', DrawParams.lw)
    plot(SimParams.time, zeros(1,SimParams.nT), 'k--')
    legend('dxL/dt', '\chi_idf/dz', '\chi_idf/dz(z_i)', '-\mu_idln(\rho_i)/dz', 'v_i')
    xlabel('time')
    ylabel('drift speed (\mum/s)')
    ylim([-5 20])
    title(['TB = ', num2str(CW.CWBias(cwInd),'%.3f')])
    
    subplot(4,1,2)
    hold on
    plot(SimParams.time, AnalysisResults.peakPhenoPos(cwInd,:) - AnalysisResults.xMaxGrad, 'r', 'LineWidth', DrawParams.lw)
    plot(SimParams.time, zeros(size(SimParams.time)), 'r--', 'LineWidth', DrawParams.lw)
    xlabel('time')
    ylabel('peak position (\mum)')
    ylim([-2 3] * 1e3)
    
    subplot(4,1,3)
    hold on
    plot(SimParams.time, AnalysisResults.rhoLinMaxGrad(cwInd,:) .* AnalysisResults.dxdtMaxGrad, 'r', 'LineWidth', DrawParams.lw)
    plot(SimParams.time, AnalysisResults.rhoLinMaxGrad(cwInd,:) .* AnalysisResults.cfDriftMaxGrad(cwInd,:), 'b', 'LineWidth', DrawParams.lw)
    plot(SimParams.time, AnalysisResults.rhoLinPhePeak(cwInd,:) .* AnalysisResults.cfDriftPhePeak(cwInd,:), 'b--', 'LineWidth', DrawParams.lw)
    plot(SimParams.time, AnalysisResults.rhoLinMaxGrad(cwInd,:) .* AnalysisResults.mrDriftMaxGrad(cwInd,:), 'm', 'LineWidth', DrawParams.lw)
    plot(SimParams.time, AnalysisResults.rhoLinMaxGrad(cwInd,:) .* AnalysisResults.relSpeedMaxGrad(cwInd,:), 'k', 'LineWidth', DrawParams.lw)
    legend('\rho_idxL/dt', '\rho_i\chi_idf/dz', '\rho_i\chi_idf/dz(z_i)', '-\mu_id\rho_i/dz', '\rho_iv_i')
    xlabel('time')
    ylabel('number flux (/s)')
    ylim([-0.1 1] * 1e2)
    
    subplot(4,1,4)
    hold on
    plot(SimParams.time, AnalysisResults.cfDivMaxGrad(cwInd,:), 'b', 'LineWidth', DrawParams.lw)
    plot(SimParams.time, AnalysisResults.mrDivMaxGrad(cwInd,:), 'm', 'LineWidth', DrawParams.lw)
    plot(SimParams.time, AnalysisResults.cfDivMaxGrad(cwInd,:) + AnalysisResults.mrDivMaxGrad(cwInd,:), 'k', 'LineWidth', DrawParams.lw)
    legend('\chi_id^2f/dz^2', '-\mu_id^2\rho_i/dz^2', 'dv_i/dz', 'location', 'SouthEast')
    xlabel('time')
    ylabel('velocity divergence (/s)')
    ylim([-8 3] * 1e-3)
    
    set(gcf, 'PaperPosition', [0 0 2 4] * DrawParams.figW);
    figName = [saveDir, 'grad_max_drift_TB=', num2str(CW.CWBias(cwInd),'%.3f'), '.png'];
    disp(['Saving figure ', figName, '...'])
    print(gcf, figName, '-dpng', ['-r',num2str(DrawParams.rez)], '-opengl')
    close(gcf)
end

end
function drawzerothordersolution(tInd, cwInd, Ui, WaveResults, SimParams, CW, MFTmodel, DrawParams, Names)
% This function draws the 0th order solution according to wave ansatz of a
% certain phenotype at a certain time.
disp('Drawing zeroth order solution ...')
xL = [-2200 1500];
CWi = CW.CWBias(cwInd);
tit = ['t = ', num2str(tInd*SimParams.outDt),' s, TB = ', num2str(CWi,'%.2f')];

figure('visible', 'off')
subplot(5,1,1)
plot(WaveResults.z, WaveResults.f)
xlabel('z (\mum)')
ylabel('perceived signal f')
xlim(xL)
title(tit)

subplot(5,1,2)
plot(WaveResults.z, WaveResults.gradient)
xlabel('z (\mum)')
ylabel('gradient (\mum^-^1)')
xlim(xL)

subplot(5,1,3)
plot(WaveResults.z, WaveResults.dynDen(cwInd,:))
hold on
[uMax, iMax] = max(WaveResults.dynDen(cwInd,:));
plot(WaveResults.z, uMax*exp((Ui(iMax)-Ui)/MFTmodel.mu(cwInd)))
xlabel('z (\mum)')
ylabel('\rho(z,TB)')
xlim(xL)
ylim([0 0.06])
legend('Simulation','Theory')

subplot(5,1,4)
plot(WaveResults.z, Ui/MFTmodel.mu(cwInd))
xlabel('z (\mum)')
ylabel('U(z,TB)')
xlim(xL)
ylim([-100 0])

subplot(5,1,5)
hold on
drift = MFTmodel.chi(cwInd) * WaveResults.gradient;
diffusion = - MFTmodel.mu(cwInd) * [0 (WaveResults.dynDen(cwInd,3:end) - WaveResults.dynDen(cwInd,1:end-2))/2/SimParams.dx./WaveResults.dynDen(cwInd,2:end-1) 0];
plot(WaveResults.z, drift)
plot(WaveResults.z, diffusion)
plot(WaveResults.z, drift+diffusion)
xlabel('z (\mum)')
ylabel('(\mum/s)')
% legend('Drift','Diffusion','Drift+Diffusion')
xlim(xL)
ylim([-2 7])

saveDir = [Names.analysisDir, 'zeroth_order_profiles\'];
if ~exist(saveDir, 'dir')
    mkdir(saveDir)
end
figName = [saveDir, tit, '.png'];
disp(['Saving figure ', figName, '...'])
set(gcf, 'PaperPosition', [0 0 0.8 1.6] * DrawParams.figW);
print(gcf, figName, '-dpng', ['-r',num2str(DrawParams.rez)], '-opengl')
close(gcf)
end
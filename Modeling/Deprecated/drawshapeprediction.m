function drawshapeprediction(tInd, WaveResults, SimParams, CW, MFTmodel, DrawParams, Names)
% This function predicts wave shape from phenotypic distribution and
% compares with simulation.

% OUTPUT
% 'predicted_profiles_t_',num2str(tInd * SimParams.outDt),'.png'

disp(['Drawing steady-state wave ansatz shape prediction ', num2str(tInd), ' ...'])
rho = sum(WaveResults.dynDen);
waveBound = DrawParams.waveBound;

% predict
WavePrediction = shapeprediction(WaveResults, SimParams, CW, MFTmodel);
z = WavePrediction.z;

% find z position of different definitions of chi min
[~, indTBAtChiMin]        = min(abs(WaveResults.chiMin - MFTmodel.chi));
zAtChiMin                 = WaveResults.zMax(indTBAtChiMin);
[~, indTBAtChiMinAtLastEps] = min(abs(WaveResults.chiMinAtLastEps - MFTmodel.chi));
zAtChiMinAtLastEps          = WaveResults.zMax(indTBAtChiMinAtLastEps);

% show predicted phenotype locations
figure('visible', 'off')
subplot(7,1,1)
hold on
plot(WaveResults.zMax, MFTmodel.chi, 'b-', 'LineWidth', DrawParams.lw)
plot(WavePrediction.zMaxPred, MFTmodel.chi, 'm-', 'LineWidth', DrawParams.lw)
plot(z, WaveResults.chiMin * ones(size(z)), 'b--', 'LineWidth', DrawParams.lw)
plot(z, WaveResults.chiMinAtArgmaxFz * ones(size(z)), 'r--', 'LineWidth', DrawParams.lw)
plot(z, WaveResults.chiMinAtLastEps * ones(size(z)), 'k--', 'LineWidth', DrawParams.lw)
plot(zAtChiMin * ones(size(MFTmodel.chi)), MFTmodel.chi, 'b--', 'LineWidth', DrawParams.lw)
plot(WaveResults.argmaxFz * ones(size(MFTmodel.chi)), MFTmodel.chi, 'r--', 'LineWidth', DrawParams.lw)
plot(zAtChiMinAtLastEps * ones(size(MFTmodel.chi)), MFTmodel.chi, 'k--', 'LineWidth', DrawParams.lw)
legend('simulatin', 'theory', ...
    'min \chi = c / max f''(z)', 'min \chi = \chi(argmax f''(z))', 'min \chi = where CDF = \epsilon', ...
    'location', 'NorthWest')
title(['t = ', num2str(tInd * SimParams.outDt), ' s'])
xlabel('peak z (\mum)')
ylabel('\chi (\mum^2/s)')
xlim(waveBound)
ylim([0 15000])

% show predicted cell density profile
subplot(7,1,2)
yRange = [0 8];
hold on
plot(z, rho, 'b', 'LineWidth', DrawParams.lw)
plot(z, WavePrediction.rhoZPred, 'm', 'LineWidth', DrawParams.lw)
plot(z, WavePrediction.rhoZGaussianSmoothedPred, 'm--', 'LineWidth', DrawParams.lw)
plot(z, WavePrediction.rhoZPotentialSmoothedPred, 'm:', 'LineWidth', DrawParams.lw)
plot(zAtChiMin * ones(size(yRange)), yRange, 'b--', 'LineWidth', DrawParams.lw)
plot(WaveResults.argmaxFz * ones(size(yRange)), yRange, 'r--', 'LineWidth', DrawParams.lw)
plot(zAtChiMinAtLastEps * ones(size(yRange)), yRange, 'k--', 'LineWidth', DrawParams.lw)
legend('simulatin', 'theory (\mu=0)', 'theory (Gaussian)', 'theory (potential)', 'location', 'NorthWest')
xlabel('z (\mum)')
ylabel('\rho (OD)')
xlim(waveBound)
ylim(yRange)

% show predicted signal
subplot(7,1,3)
yRange = [0 SimParams.asp*1.1];
hold on
plot(z, WaveResults.asp, 'b', 'LineWidth', DrawParams.lw)
plot(z, WavePrediction.aspPred, 'm', 'LineWidth', DrawParams.lw)
plot(zAtChiMin * ones(size(yRange)), yRange, 'b--', 'LineWidth', DrawParams.lw)
plot(WaveResults.argmaxFz * ones(size(yRange)), yRange, 'r--', 'LineWidth', DrawParams.lw)
plot(zAtChiMinAtLastEps * ones(size(yRange)), yRange, 'k--', 'LineWidth', DrawParams.lw)
legend('simulatin', 'theory', 'location', 'NorthWest')
xlabel('z (\mum)')
ylabel('A (\muM)')
xlim(waveBound)
ylim(yRange)

% show predicted gradient
subplot(7,1,4)
yRange = [0 3e-3];
hold on
plot(z, WaveResults.gradient, 'b', 'LineWidth', DrawParams.lw)
plot(z, WavePrediction.dfdzPred, 'm', 'LineWidth', DrawParams.lw)
plot(z, WavePrediction.dfdzPredCorrected, 'm--', 'LineWidth', DrawParams.lw)
plot(zAtChiMin * ones(size(yRange)), yRange, 'b--', 'LineWidth', DrawParams.lw)
plot(WaveResults.argmaxFz * ones(size(yRange)), yRange, 'r--', 'LineWidth', DrawParams.lw)
plot(zAtChiMinAtLastEps * ones(size(yRange)), yRange, 'k--', 'LineWidth', DrawParams.lw)
legend('simulatin', 'theory', 'theory (corrected)', 'location', 'NorthWest')
xlabel('z (\mum)')
ylabel('df/dz (/\mum)')
xlim(waveBound)
ylim(yRange)

% show predicted 2nd order gradient
subplot(7,1,5)
yRange = [-2e-6 2e-6];
hold on
plot(z, WaveResults.f2nd, 'b', 'LineWidth', DrawParams.lw)
plot(z, WavePrediction.d2fdz2Pred, 'm', 'LineWidth', DrawParams.lw)
plot(zAtChiMin * ones(size(yRange)), yRange, 'b--', 'LineWidth', DrawParams.lw)
plot(WaveResults.argmaxFz * ones(size(yRange)), yRange, 'r--', 'LineWidth', DrawParams.lw)
plot(zAtChiMinAtLastEps * ones(size(yRange)), yRange, 'k--', 'LineWidth', DrawParams.lw)
legend('simulatin', 'theory', 'location', 'NorthWest')
xlabel('z (\mum)')
ylabel('df^2/dz^2 (/\mum^2)')
xlim(waveBound)
ylim(yRange)

% show predicted speed and speed limit
subplot(7,1,6)
yRange = [-1 1];
hold on
plot(z, WaveResults.oxy, 'b', 'LineWidth', DrawParams.lw)
plot(z, WavePrediction.oxyPredExact, 'm', 'LineWidth', DrawParams.lw)
plot(z, WavePrediction.oxyPredApproximate, 'm--', 'LineWidth', DrawParams.lw)
plot(zAtChiMin * ones(size(yRange)), yRange, 'b--', 'LineWidth', DrawParams.lw)
plot(WaveResults.argmaxFz * ones(size(yRange)), yRange, 'r--', 'LineWidth', DrawParams.lw)
plot(zAtChiMinAtLastEps * ones(size(yRange)), yRange, 'k--', 'LineWidth', DrawParams.lw)
legend('simulatin', 'theory (exact)', 'theory (approximate)', 'location', 'NorthWest')
axis tight
xlabel('z (\mum)')
ylabel('O (\muM)')
xlim(waveBound)
if SimParams.oxy > 0
    yRange = [0 SimParams.oxy];
end
ylim(yRange)

% show predicted oxygen profile
subplot(7,1,7)
yRange = [0 20];
hold on
plot(z, WaveResults.speedLimit, 'b', 'LineWidth', DrawParams.lw/4)
plot(z, WavePrediction.speedLimitPred, 'm', 'LineWidth', DrawParams.lw/4)
plot(z, ones(size(z)) * WaveResults.c, 'b:', 'LineWidth', DrawParams.lw/2)
plot(z, ones(size(z)) * WavePrediction.cPred, 'm:', 'LineWidth', DrawParams.lw/2)
plot(WaveResults.z(WaveResults.phenotypeMaxInd), WaveResults.gradient(WaveResults.phenotypeMaxInd)'.*MFTmodel.chi, 'b--', 'LineWidth', DrawParams.lw)
plot(z, WavePrediction.chiOfZInterp .* WavePrediction.dfdzPred, 'm--', 'LineWidth', DrawParams.lw)
plot(zAtChiMin * ones(size(yRange)), yRange, 'b--', 'LineWidth', DrawParams.lw)
plot(WaveResults.argmaxFz * ones(size(yRange)), yRange, 'r--', 'LineWidth', DrawParams.lw)
plot(zAtChiMinAtLastEps * ones(size(yRange)), yRange, 'k--', 'LineWidth', DrawParams.lw)
legend('simulatin', 'theory', 'c (peak speed)', '\alpha_AN/A_0 (prediction)', '\chif'' (simulation)', '\chif'' (theory)', 'location', 'NorthWest')
axis tight
xlabel('z (\mum)')
ylabel('speed limit (\mum/s)')
xlim(waveBound)
ylim(yRange)

set(gcf, 'PaperPosition', [0 0 2 5] * DrawParams.figW);
print(gcf, [Names.analysisDir, 'predicted_profiles_t_',num2str(tInd * SimParams.outDt),'.png'], '-dpng', ['-r',num2str(DrawParams.rez)], '-opengl')
close(gcf)
end
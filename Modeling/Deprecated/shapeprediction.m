function WavePrediction = shapeprediction(WaveResults, SimParams, CW, MFTmodel)
% This function predicts wave shape from phenotypic distribution.

% OUTPUT_FIELDS
% z                            moving coordinate in um
% cPred                        predicted wave speed using consumption rate, cell number and total concentration
% zMaxPred                     predicted peak position in the moving coordinate z for each phenotype in um
% boundaryZInd                 index in z indicating where the maximal gradient is predicted to be (at where CDF = epsilon)
% chiOfZInterp                 predicted chemotactic coefficient in um^2/s of the phenotype peaked at z
% rhoZPred                     predicted cell density distribution in OD assuming mu = 0
% rhoZGaussianSmoothedPred     predicted cell density distribution in OD smoothed by Gaussians
% rhoZPotentialSmoothedPred    predicted cell density distribution in OD smoothed by estimated potentials
% aspPred                      predicted aspartate concentration in uM
% fPred                        predicted perceived signal
% dfdzPred                     predicted perceived gradient in um^-1
% d2fdz2Pred                   predicted perceived gradient second derivative in um^-2
% oxyPredExact                 predicted oxygen concentration in uM using the exact exponential kernel solution
% oxyPredApproximate           predicted oxygen concentration in uM using the approximate solution assuming fast consumption
% speedLimitPred               predicted upper bound of the traveling speed in um/s in the wave ansatz to ensure a phenotype locally has a peak (wave_speed_dispersion.nb)
% slowdownZInd                 index in z indicating where the wave starts to slow down, which is where gradient is predicted to be flattend
% dfdzPredCorrected            predicted perceived gradient in um^-1, after correcting leakage by making the gradient flat beyond slowdownZInd

% initialization
% disp('Calculating wave prediction ...')
z         = WaveResults.z;
% c         = WaveResults.c;
tbDistAll = WaveResults.tbDistAll;
Ncells    = sum(sum(WaveResults.dynDen))*SimParams.dx;
cPred     = Ncells * SimParams.aA / SimParams.asp;
WavePrediction = [];
if isempty(WaveResults.dynDen)
    disp(['No wave at frame ', num2str(tInd)])
    return
end

% predict zMax from tbDist
tbCDFReverse = cumsum(tbDistAll, 'reverse')*CW.dCWBias;
integrand = tbDistAll.*MFTmodel.chi./(tbCDFReverse+SimParams.KiA/SimParams.asp); % effect_of_Ki_and_P(chi).nb
integrand(isnan(integrand)) = 0;
zMaxPred = cumsum(integrand, 'reverse')*CW.dCWBias/cPred; % wave_formation_mu=0_DA=0_KA=0.nb
% includedInds = find(CW.CWBias < WaveResults.tbMaxAtChiMin); % only consider phenotypes within the wave to adjust z0
includedInds = find(CW.CWBias < WaveResults.tbMaxAtLastEps); % only consider phenotypes within the wave to adjust z0
z0 = 0;
if isfield(WaveResults, 'zMax') && ~isempty(includedInds)
    z0 = mean(WaveResults.zMax(includedInds) - zMaxPred(includedInds)); % determine integration constant to compare with simulation
end
zMaxPred = zMaxPred + z0;
if ~isempty(includedInds)
    boundaryTBInd = includedInds(end);
else
    boundaryTBInd = 1;
end
[~, boundaryZInd] = min(abs(zMaxPred(boundaryTBInd) - z));

% use predicted z of chi to find chi as a function of z, and hence aspartate profile
interpInds = find((diff([zMaxPred; zMaxPred(end)]) ~= 0) & ~isnan(zMaxPred)); % by construction zMaxPred is non-increasing, we remove where it's not changing
chiOfZInterp = interp1(zMaxPred(interpInds), MFTmodel.chi(interpInds), z);
CDFOfZInterp = interp1(zMaxPred(interpInds), tbCDFReverse(interpInds), z);
dfdzPred = cPred./chiOfZInterp;
dfdzPred(isnan(dfdzPred)) = 0;
fPred = cumsum(dfdzPred)*SimParams.dx;
% efPred = exp(fPred);
% aspPred = (SimParams.asp+SimParams.KiA)/max(efPred)*efPred-SimParams.KiA;
% aspPred = SimParams.asp/(max(efPred) - min(efPred))*efPred-SimParams.KiA;
% aspPred = (exp(fPred)-1)*SimParams.KiA;
% aspPred = efPred/max(efPred)*SimParams.asp;
aspPred = SimParams.asp * CDFOfZInterp;
aspPred = fillmissing(aspPred, 'nearest');


% use predicted z of chi to find the cumulative distribution of phenotypes as a function of z
rhoZPred = dfdzPred .* (CDFOfZInterp+SimParams.KiA/SimParams.asp) * Ncells;
rhoZPred(isnan(rhoZPred)) = 0;

% iteratively redefine gradient for better prediction
% dfdzPred = SimParams.aA/cPred * aspPred./(aspPred+SimParams.KA)./(aspPred+SimParams.KiA) .* rhoZPred; % failed

% use Gaussian approximation to find spatial density distribution when mu > 0
rhoZGaussianSmoothedPred = zeros(size(rhoZPred));
muOfZInterp = interp1(zMaxPred(interpInds), MFTmodel.mu(interpInds), z);
d2fdz2Pred = [0 (dfdzPred(3:end) - dfdzPred(1:end-2)) / SimParams.dx / 2 0];
% figure
for i = numel(chiOfZInterp):-1:1
    sigma2 = SimParams.dx^2;
    if ~(isnan(rhoZPred(i)) || rhoZPred(i) == 0 || isnan(muOfZInterp(i)) || isnan(chiOfZInterp(i)) || isnan(d2fdz2Pred(i)) || d2fdz2Pred(i) >= 0)
%         rhoZGaussianSmoothedPred(i) = rhoZGaussianSmoothedPred(i) + rhoZPred(i);
%         continue
        sigma2 = - muOfZInterp(i)/chiOfZInterp(i)/d2fdz2Pred(i); % backward update sigma2 and keep its value for last few phenotypes
    end
    gKernel = exp(- (z - z(i)).^2 / (2*sigma2)) / sqrt(2*pi*sigma2) * SimParams.dx;
    gKernel(isnan(gKernel)) = 0;
    rhoZGaussianSmoothedPred = rhoZGaussianSmoothedPred + gKernel * rhoZPred(i);
%     subplot(2,1,1)
%     hold on
%     plot(gKernel * rhoZPred(i))
end

% use potential solution to find spatial density distribution when mu > 0
rhoZPotentialSmoothedPred = zeros(size(rhoZPred));
muOfZInterp = interp1(zMaxPred(interpInds), MFTmodel.mu(interpInds), z);
% d2fdz2Pred = [0 (dfdzPred(3:end) - dfdzPred(1:end-2)) / SimParams.dx / 2 0];
for i = numel(chiOfZInterp):-1:1
    sigma2 = SimParams.dx^2;
    pKernel = exp(- (z - z(end)).^2 / (2*sigma2)) / sqrt(2*pi*sigma2) * SimParams.dx;
    if ~(isnan(rhoZPred(i)) || rhoZPred(i) == 0 || isnan(muOfZInterp(i)) || isnan(chiOfZInterp(i)) || isnan(d2fdz2Pred(i)) || d2fdz2Pred(i) >= 0)
%         Ui = (c*z - chiOfZInterp(i)*WaveResults.f);
        Ui = (cPred*z - chiOfZInterp(i)*fPred);
        [pks,locs] = findpeaks(Ui);
        if ~isempty(locs)
            Ui(1:locs) = pks;
            pKernel = exp(-Ui/muOfZInterp(i));
            pKernel(isnan(pKernel)) = 0; % backward update the potential kernel and only shift it for last few phenotypes
        end
    else
        pKernel = [pKernel(2:end) 0];
    end
    pKernel = pKernel/sum(pKernel); 
    rhoZPotentialSmoothedPred = rhoZPotentialSmoothedPred + pKernel * rhoZPred(i);
%     subplot(2,1,2)
%     hold on
%     plot(pKernel * rhoZPred(i))
end

% use predicted cell density to find oxygen profile
if SimParams.oxy > 0
    oxyPredExact = SimParams.oxy - SimParams.aO/cPred * exp(SimParams.kO/cPred * z) .* cumsum(exp(- SimParams.kO/cPred * z) .* rhoZPred, 'reverse') * SimParams.dx;
    oxyPredApproximate = SimParams.oxy - SimParams.aO/SimParams.kO * rhoZPred;
    oxyPredExact(oxyPredExact < 0) = 0;
    oxyPredApproximate(oxyPredApproximate < 0) = 0;
else
    oxyPredExact = zeros(size(rhoZPred));
    oxyPredApproximate = zeros(size(rhoZPred));
end

% find speed limit and use it to predict where wave starts to slow down
speedLimitPred = SimParams.aA * (aspPred.^2 - SimParams.KiA * SimParams.KA) ...
    ./(aspPred + SimParams.KA).^2./(aspPred + SimParams.KiA) ...
    .* rhoZPred.^2 ...
    ./ [0 (rhoZPred(3:end) - rhoZPred(1:end-2))/(2*SimParams.dx) 0];
indSpeedLimit = find(speedLimitPred > cPred);
if ~isempty(indSpeedLimit)
    slowdownZInd = indSpeedLimit(1);
else
    slowdownZInd = 1;
end
dfdzPredCorrected = dfdzPred;
dfdzPredCorrected(1:slowdownZInd) = dfdzPredCorrected(slowdownZInd);

% explore prediction discrepancy
% subplot(4,1,1)
% hold on
% plot(z, speedLimitPred)
% plot(z, cPred*ones(size(z)))
% ylim([0 15])
% subplot(4,1,2)
% hold on
% plot(z, dfdzPred)
% plot(z, SimParams.aA/cPred * aspPred./(aspPred+SimParams.KA)./(aspPred+SimParams.KiA) .* rhoZPred)
% subplot(4,1,3)
% hold on
% plot(z, [0 (aspPred(3:end) - aspPred(1:end-2))/(2*SimParams.dx) 0])
% plot(z, SimParams.aA/cPred * aspPred./(aspPred+SimParams.KA) .* rhoZPred)
% subplot(4,1,4)
% hold on
% plot(z, aspPred)
% plot(z, SimParams.aA/cPred * cumsum(rhoZPred) * SimParams.dx)

% save prediction
WavePrediction.z                         = z;
WavePrediction.cPred                     = cPred;
WavePrediction.zMaxPred                  = zMaxPred;
WavePrediction.boundaryZInd              = boundaryZInd;
WavePrediction.chiOfZInterp              = chiOfZInterp;
WavePrediction.rhoZPred                  = rhoZPred;
WavePrediction.rhoZGaussianSmoothedPred  = rhoZGaussianSmoothedPred;
WavePrediction.rhoZPotentialSmoothedPred = rhoZPotentialSmoothedPred;
WavePrediction.aspPred                   = aspPred;
WavePrediction.fPred                     = fPred;
WavePrediction.dfdzPred                  = dfdzPred;
WavePrediction.d2fdz2Pred                = d2fdz2Pred;
WavePrediction.oxyPredExact              = oxyPredExact;
WavePrediction.oxyPredApproximate        = oxyPredApproximate;
WavePrediction.speedLimitPred            = speedLimitPred;
WavePrediction.slowdownZInd              = slowdownZInd;
WavePrediction.dfdzPredCorrected         = dfdzPredCorrected;
end
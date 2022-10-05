function WaveResults = analyzewavesnapshot(tInd, wb, wf, AnalysisParams, SimResults, SimParams, CW, MFTmodel)
% This function analyzes wave peak properties in total and for each phenotype
% at a specified time.

% INPUTS
% tInd : time index to analyze wave
% wb   : index length from the back to the peak of the wave at this time
% wf   : index length from the front to the peak of the wave at this time

% OUTPUT_FIELDS
% z                   moving coordinate in um
% dynDen              CW.CWBias * z array of cell density in OD in moving coordinate z at time t, starting from minimal total cell density
% asp                 asp concentration in uM in moving coordinate z at time t
% oxy                 oxy concentration in uM in moving coordinate z at time t
% glu                 glu concentration in uM in moving coordinate z at time t
% ser                 ser concentration in uM in moving coordinate z at time t
% nutrient            nutrient concentration in uM in moving coordinate z at time t
% f                   perceived signal in moving coordinate z at time t
% fMax                maximum of perceived signal (used for normalization)
% gradient            perceived gradient in 1/um in moving coordinate z at time t
% gradientMax         maximum of perceived gradient in 1/um
% gradientMaxInd      z index of the maximum of perceived gradient
% xGradientMax        x position of the maximum of perceived gradient in um
% rhoLGradientMax     phenotype linear density in cells/um at the maximum of perceived gradient
% chemGradientMax     phenotype chemotactic drift term in um/s at the maximum of perceived gradient
% diffGradientMax     phenotype diffusive drift term in um/s at the maximum of perceived gradient
% chemDivGradMax      divergence of phenotype chemotactic drift term in 1/s at the maximum of perceived gradient
% diffDivGradMax      divergence of phenotype diffusive drift term in 1/s at the maximum of perceived gradient
% f2nd                perceived gradient change in 1/um^2 in moving coordinate z at time t
% f2ndMax             maximum of perceived gradient change in 1/um^2
% argmaxFz            z in um where gradientMax is
% chiMin              the minimal chi in um^2/s in the wave from c = chi*gradientMax
% zMax                a column of peak z in um for each phenotype
% tbMaxAtArgmaxFz     the maximal tb in the wave using gradientMax
% chiMinAtArgmaxFz    the minimal chi in the wave using gradientMax
% muMinAtArgmaxFz     the minimal mu in the wave using gradientMax
% tbMaxAtChiMin       the maximal tb in the wave using chiMin
% muMinAtChiMin       the minimal mu in the wave using chiMin
% tbMaxAtLastEps      the maximal tb in the wave using CDF = epsilon
% chiMinAtLastEps     the minimal chi in the wave using CDF = epsilon
% muMinAtLastEps      the minimal mu in the wave using CDF = epsilon
% tbMaxAtG0           the maximal tb in the wave using chi peaked at g0 = c / predicted f' at chi(CDF = epsilon)
% chiMinAtG0          the minimal chi in the wave using chi peaked at g0 = c / predicted f' at chi(CDF = epsilon)
% muMinAtG0           the minimal mu in the wave using chi peaked at g0 = c / predicted f' at chi(CDF = epsilon)
% chiMinAtAKiA        the minimal chi in the wave using A >> SimParams.KiA
% tbMaxAtAKiA         the maximal tb in the wave using A >> SimParams.KiA
% chiMinAtAKiAKA      the minimal chi in the wave using A >> sqrt(SimParams.KiA * SimParams.KA)
% tbMaxAtAKiAKA       the maximal tb in the wave using A >> sqrt(SimParams.KiA * SimParams.KA)
% chiMinAtF2          the minimal chi in the wave using |d^2f/dz^2| << (df/dz)^2
% tbMaxAtF2           the maximal tb in the wave using |d^2f/dz^2| << (df/dz)^2
% c                   wave speed in um/s
% denMax              peak density in OD
% denMaxInd           x index of where peak is
% phenotypeMax        phenotype peak density in OD
% phenotypeMaxInd     z index of where phenotype peak is
% phenotypeMin        phenotype back min density in OD
% xPhenotypeMax       x position in um of the phenotype peak density
% gradPhenotypeMax    perceived gradient in 1/um at the phenotype peak density
% rhoLPhenotypeMax    phenotype linear density in cells/um at the phenotype peak density
% chemPhenotypeMax    phenotype chemotactic drift term in um/s at the phenotype peak density
% waveBack            x index of where z = 1
% tbDistGMax          tumble bias distribution in the wave (defined as starting from gradient maximum)
% tbDistAll           tumble bias distribution in the wave (defined as starting from density minimum)
% speedLimit          upper bound of the traveling speed in um/s in the wave ansatz to ensure a phenotype locally has a peak (wave_speed_dispersion.nb)

% initialize results
disp(['Analyzing wave snapshot ', num2str(tInd), ' ...'])
WaveResults = [];
WaveResults.tInd = tInd;
WaveResults.wb = wb;
WaveResults.wf = wf;

% determine peak
cellDensity = squeeze(sum(SimResults.dynProfile.cellDensity(:, :, tInd))); % cell density in OD totaling all phenotypes
maxInd = numel(cellDensity);
%     wave_offset = max(1, wo);
%     [denMax, denMaxInd] = max(cell_density(wave_offset:end));
%     if isempty(denMaxInd) % skips when wave_offset is beyond channel end
%         return
%     end
%     waveBack    = max(1, wave_offset + denMaxInd - wb);
%     wave_front   = min(max_ind, wave_offset + denMaxInd + wf);
%     z = x - x(denMaxInd + wave_offset - waveBack);
dsd = diff(sign( diff(cellDensity) ));
denMaxInd = find(dsd == -2) + 1; % find where first derivative = 0 and second derivative < 0
if isempty(denMaxInd)
    denMaxInd = 1;
end
denMaxInd = denMaxInd(end); % the last local maximum is real
denMax = cellDensity(denMaxInd);

% determine wave packet range
denMinInd = find(dsd == 2) + 1; % find where first derivative = 0 and second derivative > 0
if isempty(denMinInd)
    denMinInd = 1;
end
denMinInd = denMinInd(1); % the first local minimum is real
waveBack  = max(denMinInd, denMaxInd - wb);
waveFront = min(maxInd, denMaxInd + wf);

% calculate moving coordinate
dynDen = squeeze(SimResults.dynProfile.cellDensity(:, waveBack:waveFront, tInd));
if isempty(dynDen)
    WaveResults.dynDen = dynDen;
    disp(['No wave at frame ', num2str(tInd)])
    return
end
[phenotypeMax, phenotypeMaxInd] = max(dynDen, [], 2);
z1 = 0:SimParams.dx:(length(waveBack:waveFront-1)*SimParams.dx);
z = z1 - z1(max(denMaxInd - waveBack + 1, 1));

% determine each phenotype's back density
phenotypeMin = zeros(size(phenotypeMax));
for cwInd = 1:CW.ncwb
    dsd = diff(sign( diff(SimResults.dynProfile.cellDensity(cwInd, :, tInd)) ));
    denMinInd = find(dsd == 2) + 1; % find where first derivative = 0 and second derivative > 0
    if isempty(denMinInd)
        denMinInd = 1;
    end
    denMinInd = denMinInd(1); % the first local minimum is real
    phenotypeMin(cwInd) = SimResults.dynProfile.cellDensity(cwInd, denMinInd, tInd);
end

%     plot(sum(SimResults.dynProfile.cellDensity(:, :, t)))
%     plot(sum(dynDen))
%     title(['t = ', num2str(t*SimParams.outDt/60), ' min'])
%     pause

% check aspartate and O2 profile
%     KiO1 = 40;
%     KiO2 = 150;
%     KaO1 = 330;
%     KaO2 = 10000;
asp = SimResults.dynProfile.asp(waveBack:waveFront,tInd)';
%     O_profile   = dyn_profile.O(waveBack:wave_front, t);
f = log((1+asp/SimParams.KiA)./(1+asp/SimParams.KaA));
fMax = max(f);
%     fO = (log((1+O_profile/KiO1)./(1+O_profile/KaO1))-log((1+O_profile/KiO2)./(1+O_profile/KaO2)))/5000;

% find the gradient
gradient = [0 (f(3:end)-f(1:end-2))/2/SimParams.dx 0]; % gradientO=((f(3:end)-f(1:end-2))/2/SimParams.dx+(fO(3:end)-fO(1:end-2))/2/SimParams.dx);
[gradientMax, gradientMaxInd] = max(gradient);
absGradientMaxInd = waveBack-1+gradientMaxInd;
xGradientMax = SimParams.x(absGradientMaxInd);

% calculate the wave phenotype distribution
tbDistGMax = sum(dynDen(:,gradientMaxInd:end),2);
tbDistGMax = tbDistGMax/sum(tbDistGMax)/CW.dCWBias;
tbDistAll  = sum(dynDen,2);
tbDistAll  = tbDistAll/sum(tbDistAll)/CW.dCWBias;

% find each phenotype's drifts at the gradient maximum
rhoLGradientMax = dynDen(:,gradientMaxInd) * SimParams.OD2Cell / SimParams.dx;
chemGradientMax = MFTmodel.chi * gradientMax;
if absGradientMaxInd == 1
    dlnRhodz = SimResults.dynProfile.cellDensity(:, absGradientMaxInd+1, tInd) - SimResults.dynProfile.cellDensity(:, absGradientMaxInd, tInd);
elseif absGradientMaxInd == SimParams.L
    dlnRhodz = SimResults.dynProfile.cellDensity(:, absGradientMaxInd, tInd) - SimResults.dynProfile.cellDensity(:, absGradientMaxInd-1, tInd);
else
    dlnRhodz = (SimResults.dynProfile.cellDensity(:, absGradientMaxInd+1, tInd) - SimResults.dynProfile.cellDensity(:, absGradientMaxInd-1, tInd))/2;
end
dlnRhodz = squeeze(dlnRhodz ./ SimResults.dynProfile.cellDensity(:, absGradientMaxInd, tInd) / SimParams.dx);
diffGradientMax = -MFTmodel.mu .* dlnRhodz;

% find f''(z)
f2nd = [0 (f(3:end)+f(1:end-2)-2*f(2:end-1))/SimParams.dx^2 0];
f2ndMax = max(abs(f2nd));

% find each phenotype's speed divergences at the gradient maximum
chemDivGradMax = MFTmodel.chi * f2nd(gradientMaxInd);
d2lnRhodz2 = [zeros(CW.ncwb,1) (log(dynDen(:,3:end))+log(dynDen(:,1:end-2))-2*log(dynDen(:,2:end-1)))/SimParams.dx^2 zeros(CW.ncwb,1)];
diffDivGradMax = -MFTmodel.mu .* d2lnRhodz2(:, gradientMaxInd);

if tInd > AnalysisParams.calculatePeakSpacing
    % calculate peak speed
    profile1 = sum(SimResults.dynProfile.cellDensity(:,waveBack:waveFront,tInd));
    profile2 = sum(SimResults.dynProfile.cellDensity(:,waveBack:waveFront,tInd-AnalysisParams.calculatePeakSpacing));
    [~, ind1] = max(profile1);
    [~, ind2] = max(profile2);
    c = (ind1-ind2)*SimParams.dx/AnalysisParams.calculatePeakSpacing/SimParams.outDt;
else
    c = nan;
end

% find each phenotype's peak and the boundary phenotype
zMax = z(phenotypeMaxInd)'; % a column of peak z for each phenotype
argmaxFz = z(gradientMaxInd); % z where gradientMax is
chiMin = c / gradientMax; % the minimal chi in the wave from c = chi*gradientMax

[~, indTBMaxAtArgmaxFz] = min(abs(zMax - argmaxFz));
tbMaxAtArgmaxFz = CW.CWBias(indTBMaxAtArgmaxFz); % the maximal tb in the wave using gradientMax
chiMinAtArgmaxFz = MFTmodel.chi(indTBMaxAtArgmaxFz); % the minimal chi in the wave using gradientMax
muMinAtArgmaxFz = MFTmodel.mu(indTBMaxAtArgmaxFz); % the minimal mu in the wave using gradientMax

if isempty(chiMin)
    tbMaxAtChiMin = tbMaxAtArgmaxFz;
    muMinAtChiMin = muMinAtArgmaxFz;
else
    [~, indTBMaxAtChiMin] = min(abs(MFTmodel.chi - chiMin));
    tbMaxAtChiMin = CW.CWBias(indTBMaxAtChiMin); % the maximal tb in the wave using chiMin
    muMinAtChiMin = MFTmodel.mu(indTBMaxAtChiMin); % the minimal mu in the wave using chiMin
end

tbCDF = cumsum(tbDistAll, 'reverse') * CW.dCWBias;
indInWave = find(tbCDF > SimParams.KiA/SimParams.asp);
indAtLastEpsInWave = indInWave(end);
tbMaxAtLastEps = CW.CWBias(indAtLastEpsInWave); % the maximal tb in the wave using the last phenotype in the wave
chiMinAtLastEps = MFTmodel.chi(indAtLastEpsInWave); % the minimal chi in the wave using the last phenotype in the wave
muMinAtLastEps = MFTmodel.mu(indAtLastEpsInWave); % the minimal mu in the wave using the last phenotype in the wave

% find each phenotype's drifts at the gradient maximum
xPhenotypeMax    = SimParams.x(waveBack - 1 + phenotypeMaxInd)';
gradPhenotypeMax = gradient(phenotypeMaxInd)';
rhoLPhenotypeMax = diag(dynDen(:,phenotypeMaxInd)) * SimParams.OD2Cell / SimParams.dx;
chemPhenotypeMax = MFTmodel.chi .* gradPhenotypeMax;

% find speed limit
rho = sum(dynDen);
speedLimit = SimParams.aA * (asp.^2 - SimParams.KiA * SimParams.KA) ...
    ./(asp + SimParams.KA).^2./(asp + SimParams.KiA) ...
    .* rho.^2 ...
    ./ [0 (rho(3:end) - rho(1:end-2))/(2*SimParams.dx) 0];

% save results
WaveResults.z                = z;
WaveResults.dynDen           = dynDen;
WaveResults.asp              = asp;
WaveResults.oxy              = SimResults.dynProfile.oxy(waveBack:waveFront,tInd)';
WaveResults.glu              = SimResults.dynProfile.glu(waveBack:waveFront,tInd)';
WaveResults.ser              = SimResults.dynProfile.ser(waveBack:waveFront,tInd)';
WaveResults.nutrient         = SimResults.dynProfile.oxy(waveBack:waveFront,tInd)';
WaveResults.f                = f;
WaveResults.fMax             = fMax;
WaveResults.gradient         = gradient;
WaveResults.gradientMax      = gradientMax;
WaveResults.gradientMaxInd   = gradientMaxInd;
WaveResults.xGradientMax     = xGradientMax;
WaveResults.rhoLGradientMax  = rhoLGradientMax;
WaveResults.chemGradientMax  = chemGradientMax;
WaveResults.diffGradientMax  = diffGradientMax;
WaveResults.chemDivGradMax   = chemDivGradMax;
WaveResults.diffDivGradMax   = diffDivGradMax;
WaveResults.f2nd             = f2nd;
WaveResults.f2ndMax          = f2ndMax;
WaveResults.argmaxFz         = argmaxFz;
WaveResults.chiMin           = chiMin;
WaveResults.zMax             = zMax;
WaveResults.tbMaxAtArgmaxFz  = tbMaxAtArgmaxFz;
WaveResults.chiMinAtArgmaxFz = chiMinAtArgmaxFz;
WaveResults.muMinAtArgmaxFz  = muMinAtArgmaxFz;
WaveResults.tbMaxAtChiMin    = tbMaxAtChiMin;
WaveResults.muMinAtChiMin    = muMinAtChiMin;
WaveResults.tbMaxAtLastEps   = tbMaxAtLastEps;
WaveResults.chiMinAtLastEps  = chiMinAtLastEps;
WaveResults.muMinAtLastEps   = muMinAtLastEps;
WaveResults.c                = c;
WaveResults.denMax           = denMax;
WaveResults.denMaxInd        = denMaxInd;
WaveResults.phenotypeMax     = phenotypeMax;
WaveResults.phenotypeMaxInd  = phenotypeMaxInd;
WaveResults.phenotypeMin     = phenotypeMin;
WaveResults.xPhenotypeMax    = xPhenotypeMax;
WaveResults.gradPhenotypeMax = gradPhenotypeMax;
WaveResults.rhoLPhenotypeMax = rhoLPhenotypeMax;
WaveResults.chemPhenotypeMax = chemPhenotypeMax;
WaveResults.waveBack         = waveBack;
WaveResults.tbDistGMax       = tbDistGMax;
WaveResults.tbDistAll        = tbDistAll;
WaveResults.speedLimit       = speedLimit;

% use theoretical prediction to find the boundary phenotype
if (isnan(c) || ~c)
    disp(['Warning: invalid wave speed for shape prediction calculations in analyzewavesnapshot at time = ', num2str(tInd * SimParams.outDt)])
    indTBMaxAtAKiA = 1;
    indTBMaxAtAKiAKA = 1;
    indTBMaxAtF2 = 1;
else
    WavePrediction = shapeprediction(WaveResults, SimParams, CW, MFTmodel);
    
    % use prediction to find boundary chi
    [~, indZAKiA] = min(abs(WavePrediction.aspPred - 20 * SimParams.KiA));
    zAKiA = z(indZAKiA); % z where A >> SimParams.KiA
    [~, indTBMaxAtAKiA] = min(abs(zMax - zAKiA));
    
    [~, indZAKiAKA] = min(abs(WavePrediction.aspPred - 100 * sqrt(SimParams.KiA * SimParams.KA)));
    zAKiAKA = z(indZAKiAKA); % z where A >> sqrt(SimParams.KiA * SimParams.KA)
    [~, indTBMaxAtAKiAKA] = min(abs(zMax - zAKiAKA));
    
    [~, indZF2] = min(abs(-WavePrediction.d2fdz2Pred./WavePrediction.dfdzPred.^2 - 0.25));
    zF2 = z(indZF2); % z where |d^2f/dz^2| << (df/dz)^2
    [~, indTBMaxAtF2] = min(abs(zMax - zF2));
    
%     g0 = WavePrediction.cPred / chiMinAtLastEps;
end
tbMaxAtAKiA = CW.CWBias(indTBMaxAtAKiA); % the maximal tb in the wave using A >> SimParams.KiA
chiMinAtAKiA = MFTmodel.chi(indTBMaxAtAKiA); % the minimal chi in the wave using A >> SimParams.KiA
tbMaxAtAKiAKA = CW.CWBias(indTBMaxAtAKiAKA); % the maximal tb in the wave using A >> sqrt(SimParams.KiA * SimParams.KA)
chiMinAtAKiAKA = MFTmodel.chi(indTBMaxAtAKiAKA); % the minimal chi in the wave using A >> sqrt(SimParams.KiA * SimParams.KA)
tbMaxAtF2 = CW.CWBias(indTBMaxAtF2); % the maximal tb in the wave using |d^2f/dz^2| << (df/dz)^2
chiMinAtF2 = MFTmodel.chi(indTBMaxAtF2); % the minimal chi in the wave using |d^2f/dz^2| << (df/dz)^2

WaveResults.chiMinAtAKiA     = chiMinAtAKiA;
WaveResults.tbMaxAtAKiA      = tbMaxAtAKiA;
WaveResults.chiMinAtAKiAKA   = chiMinAtAKiAKA;
WaveResults.tbMaxAtAKiAKA    = tbMaxAtAKiAKA;
WaveResults.chiMinAtF2       = chiMinAtF2;
WaveResults.tbMaxAtF2        = tbMaxAtF2;
end
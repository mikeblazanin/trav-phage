function AnalysisResults = analyzewavedynamics(AnalysisParams, SimResults, ...
    SimParams, CW, MFTmodel, DrawParams, Names, RunParams)
% This function analyze wave results over all times.

%% initialization
disp('Analyzing wave dynamics ...')
AnalysisResults = [];
AnalysisResults.peakPos           = nan(1,SimParams.nT); % um, position of wave peak
AnalysisResults.peakSpd           = nan(1,SimParams.nT); % um/s, speed of the wave peak
AnalysisResults.peakDen           = nan(1,SimParams.nT); % cells/ml, volume density at cell peak
AnalysisResults.waveCellsGMax     = nan(1,SimParams.nT); % number of cells in the wave (defined as starting from the maximal gradient)
AnalysisResults.waveCellsAll      = nan(1,SimParams.nT); % number of cells in the wave (defined as starting from the minimal of total cell density), this is mostly only used when predicting wave speed
AnalysisResults.gradientMax       = nan(1,SimParams.nT); % the maximal gradient in 1/um
AnalysisResults.chiMin            = nan(1,SimParams.nT); % the minimal chi in um^2/s in the wave from c = chi*gradientMax
AnalysisResults.tbMaxAtChiMin     = nan(1,SimParams.nT); % the maximal tb in the wave using chiMin
AnalysisResults.muMinAtChiMin     = nan(1,SimParams.nT); % the minimal chi in um^2/s in the wave using chiMin
AnalysisResults.tbMaxAtArgmaxFz   = nan(1,SimParams.nT); % the maximal tb in the wave using gradientMax
AnalysisResults.chiMinAtArgmaxFz  = nan(1,SimParams.nT); % the minimal chi in um^2/s in the wave using gradientMax
AnalysisResults.muMinAtArgmaxFz   = nan(1,SimParams.nT); % the minimal mu in um^2/s in the wave using gradientMax
AnalysisResults.tbMaxAtLastEps    = nan(1,SimParams.nT); % the maximal tb in the wave using CDF = epsilon
AnalysisResults.chiMinAtLastEps   = nan(1,SimParams.nT); % the minimal chi in um^2/s in the wave using CDF = epsilon
AnalysisResults.muMinAtLastEps    = nan(1,SimParams.nT); % the minimal mu in um^2/s in the wave using CDF = epsilon
AnalysisResults.chiMinPred        = nan(1,SimParams.nT); % the minimal chi in um^2/s in the wave from c = chi*max{dfdzPred}
AnalysisResults.tbMaxAtChiMinPred = nan(1,SimParams.nT); % the maximal tb in the wave using chiMinPred
AnalysisResults.muMinAtChiMinPred = nan(1,SimParams.nT); % the minimal chi in um^2/s in the wave using chiMinPred
AnalysisResults.fluxBackCells     = zeros(CW.ncwb,SimParams.nT); % number flux of cells in the back of the wave (up to max gradient) of each phenotype (Risken 5.111, Kramers)
AnalysisResults.imFluxBackCells   = zeros(CW.ncwb,SimParams.nT); % improved number flux of cells in the back of the wave (up to max gradient) of each phenotype (Risken 5.112, Kramers higher orders)
AnalysisResults.acFluxBackCells   = zeros(CW.ncwb,SimParams.nT); % more accurate number flux of cells in the back of the wave (up to max gradient) of each phenotype (Risken 5.109, non-Gaussian potential)
% AnalysisResults.cpFluxBackCells   = zeros(CW.ncwb,SimParams.nT); % c * rho at back boundary, number flux of cells in the back of the wave of each phenotype
% AnalysisResults.cqFluxBackCells   = zeros(CW.ncwb,SimParams.nT); % c * rho at each phenotype's back, number flux of cells in the back of the wave of each phenotype
% AnalysisResults.fpFluxBackCells   = zeros(CW.ncwb,SimParams.nT); % chi * df/dz at peak * rho at back boundary, number flux of cells in the back of the wave of each phenotype
% AnalysisResults.fqFluxBackCells   = zeros(CW.ncwb,SimParams.nT); % chi * df/dz at peak * rho at each phenotype's back, number flux of cells in the back of the wave of each phenotype
AnalysisResults.peakPhenoPos      = nan(CW.ncwb,SimParams.nT); % um, peak position of each phenotype
AnalysisResults.peakPhenoSpd      = nan(CW.ncwb,SimParams.nT); % um/s, peak speed of each phenotype
AnalysisResults.peakPhenoDen      = nan(CW.ncwb,SimParams.nT); % cells/ml, peak density of each phenotype
AnalysisResults.tbExclude         = ones(1,SimParams.nT); % lowest tumble bias of cells that don't have a peak in the wave
AnalysisResults.PWaveGMax         = nan(CW.ncwb,SimParams.nT); % pdf of tumble bias in the wave (defined as starting from the maximal gradient)
AnalysisResults.PWaveAll          = nan(CW.ncwb,SimParams.nT); % pdf of tumble bias in the wave (defined as starting from the minimal density)
AnalysisResults.U                 = cell(1,SimParams.nT); % potential at each time, each cell has dimension CW x WaveResults.z
AnalysisResults.xMaxGrad          = zeros(1,SimParams.nT); % um, x of max gradient
AnalysisResults.dxdtMaxGrad       = nan(1,SimParams.nT); % um/s, dx/dt of max gradient
AnalysisResults.rhoLinMaxGrad     = zeros(CW.ncwb,SimParams.nT); % cells/um, linear density of cells at max gradient
AnalysisResults.cfDriftMaxGrad    = zeros(CW.ncwb,SimParams.nT); % um/s, chi * df/dx at max gradient
AnalysisResults.mrDriftMaxGrad    = zeros(CW.ncwb,SimParams.nT); % um/s, -mu * dln(rho)/dx at max gradient
AnalysisResults.cfDivMaxGrad      = zeros(CW.ncwb,SimParams.nT); % 1/s, chi * d^2f/dx^2 at max gradient
AnalysisResults.mrDivMaxGrad      = zeros(CW.ncwb,SimParams.nT); % 1/s, -mu * d^2ln(rho)/dx^2 at max gradient
AnalysisResults.relSpeedMaxGrad   = zeros(CW.ncwb,SimParams.nT); % um/s, drift speed of each phenotype in the moving cordinate defined by dxdtMaxGrad
AnalysisResults.rhoLinPhePeak     = zeros(CW.ncwb,SimParams.nT); % cells/um, linear density of cells at each phenotype's peak
AnalysisResults.cfDriftPhePeak    = zeros(CW.ncwb,SimParams.nT); % um/s, chi * df/dx at each phenotype's peak
AnalysisResults.gsFluxBackCells   = zeros(CW.ncwb,SimParams.nT,7); % group-speed estimated number flux of cells in the back of the wave (up to max gradient) of each phenotype, 5 methods to find chi_boundary

%% find wave at each time
if RunParams.reDrawZerothOrderSol
    tIndRange = DrawParams.tInds;
    cwIndRange = 200:20:400;
else
    tIndRange = 1:SimParams.nT;
    cwIndRange = 1:CW.ncwb;
end
for tInd = tIndRange
    % analyze snapshot for aggregate results
    WaveResults = analyzewavesnapshot(tInd, AnalysisParams.wb(end), AnalysisParams.wf(end), ...
        AnalysisParams, SimResults, SimParams, CW, MFTmodel);
    if RunParams.reDrawZerothOrderSol && tInd > AnalysisParams.calculatePeakSpacing
        WaveResultsPast = analyzewavesnapshot(tInd-AnalysisParams.calculatePeakSpacing, AnalysisParams.wb(end), AnalysisParams.wf(end), ...
            AnalysisParams, SimResults, SimParams, CW, MFTmodel);
        AnalysisResults.xMaxGrad(tInd-AnalysisParams.calculatePeakSpacing) = WaveResultsPast.xGradientMax;
    end

    if isempty(WaveResults.dynDen)
        disp(['Wave is destroyed at frame ', num2str(tInd)])
        break
    end
    
    % use wave prediction to find boundary
    WavePrediction = shapeprediction(WaveResults, SimParams, CW, MFTmodel);
%     chiMinPred = WavePrediction.chiOfZInterp(WavePrediction.boundaryZInd);
    chiMinPred = WavePrediction.chiOfZInterp(WavePrediction.slowdownZInd);
    if isempty(chiMinPred)
        tbMaxAtChiMinPred = tbMaxAtChiMin;
        muMinAtChiMinPred = muMinAtChiMin;
    else
        [~, indTBMaxAtChiMinPred] = min(abs(MFTmodel.chi - chiMinPred));
        tbMaxAtChiMinPred = CW.CWBias(indTBMaxAtChiMinPred); % the maximal tb in the wave using chiMinPred
        muMinAtChiMinPred = MFTmodel.mu(indTBMaxAtChiMinPred); % the minimal mu in the wave using chiMinPred
    end
    
%     disp([num2str(WaveResults.chiMin,2), ':', num2str(WaveResults.muMinAtChiMin,2), ':', num2str(WaveResults.tbMaxAtChiMin,2)])
%     disp([num2str(chiMinPred,2), ':', num2str(muMinAtChiMinPred,2), ':', num2str(tbMaxAtChiMinPred,2)])
%     disp([num2str(WaveResults.chiMinAtLastEps,2), ':', num2str(WaveResults.muMinAtLastEps,2), ':', num2str(WaveResults.tbMaxAtLastEps,2)])
%     disp([num2str(WaveResults.chiMinAtArgmaxFz,2), ':', num2str(WaveResults.muMinAtArgmaxFz,2), ':', num2str(WaveResults.tbMaxAtArgmaxFz,2)])
%     disp([num2str(WaveResults.c,2), ':', num2str(WavePrediction.cPred,2)])
    
    % save direct results
    AnalysisResults.peakPos(tInd)            = SimParams.x(WaveResults.denMaxInd);
    AnalysisResults.peakSpd(tInd)            = WaveResults.c;
    AnalysisResults.peakDen(tInd)            = WaveResults.denMax * SimParams.OD1;
    AnalysisResults.waveCellsGMax(tInd)      = sum(sum(WaveResults.dynDen(:,WaveResults.gradientMaxInd:end))) * SimParams.OD2Cell;
    AnalysisResults.waveCellsAll(tInd)       = sum(sum(WaveResults.dynDen)) * SimParams.OD2Cell;
    AnalysisResults.gradientMax              = WaveResults.gradientMax;
    AnalysisResults.chiMin(tInd)             = WaveResults.chiMin;
    AnalysisResults.tbMaxAtChiMin(tInd)      = WaveResults.tbMaxAtChiMin;
    AnalysisResults.muMinAtChiMin(tInd)      = WaveResults.muMinAtChiMin;
    AnalysisResults.tbMaxAtArgmaxFz(tInd)    = WaveResults.tbMaxAtArgmaxFz;
    AnalysisResults.chiMinAtArgmaxFz(tInd)   = WaveResults.chiMinAtArgmaxFz;
    AnalysisResults.muMinAtArgmaxFz(tInd)    = WaveResults.muMinAtArgmaxFz;
    AnalysisResults.tbMaxAtLastEps(tInd)     = WaveResults.tbMaxAtLastEps;
    AnalysisResults.chiMinAtLastEps(tInd)    = WaveResults.chiMinAtLastEps;
    AnalysisResults.muMinAtLastEps(tInd)     = WaveResults.muMinAtLastEps;
    AnalysisResults.chiMinPred(tInd)         = chiMinPred;
    AnalysisResults.tbMaxAtChiMinPred(tInd)  = tbMaxAtChiMinPred;
    AnalysisResults.muMinAtChiMinPred(tInd)  = muMinAtChiMinPred;
    AnalysisResults.PWaveGMax(:,tInd)        = WaveResults.tbDistGMax;
    AnalysisResults.PWaveAll(:,tInd)         = WaveResults.tbDistAll;
    AnalysisResults.U{tInd}                  = nan(size(WaveResults.dynDen));
    AnalysisResults.xMaxGrad(tInd)           = WaveResults.xGradientMax;
    if tInd > AnalysisParams.calculatePeakSpacing
        AnalysisResults.dxdtMaxGrad(tInd)    = (AnalysisResults.xMaxGrad(tInd) - AnalysisResults.xMaxGrad(tInd-AnalysisParams.calculatePeakSpacing))/AnalysisParams.calculatePeakSpacing/SimParams.outDt;
    end
    AnalysisResults.rhoLinMaxGrad(:,tInd)    = WaveResults.rhoLGradientMax;
    AnalysisResults.cfDriftMaxGrad(:,tInd)   = WaveResults.chemGradientMax;
    AnalysisResults.mrDriftMaxGrad(:,tInd)   = WaveResults.diffGradientMax;
    AnalysisResults.cfDivMaxGrad(:,tInd)     = WaveResults.chemDivGradMax;
    AnalysisResults.mrDivMaxGrad(:,tInd)     = WaveResults.diffDivGradMax;
    AnalysisResults.rhoLinPhePeak(:,tInd)    = WaveResults.rhoLPhenotypeMax;
    AnalysisResults.cfDriftPhePeak(:,tInd)   = WaveResults.chemPhenotypeMax;
    
    % find where phenotype falls off
    if ~isempty(WaveResults.dynDen)
        AnalysisResults.peakPhenoPos(:,tInd) = WaveResults.xPhenotypeMax;
        AnalysisResults.peakPhenoDen(:,tInd) = WaveResults.phenotypeMax;
%                 exclude = find(phenotype_max_ind == 1);
%                 if ~isempty(exclude)
%                     tb_exclude(t) = CW.CWBias(exclude(1));
%                 end
        includeFlag = WaveResults.phenotypeMaxInd ~= 1; % whether a phenotype is included in the wave
        if ~all(includeFlag) % there are phenotypes not in the wave
            includeInd = find(includeFlag);
            if isempty(includeInd)
                bInd = 1;
            else
                bInd = includeInd(end);
                if bInd < CW.ncwb
                    bInd = bInd + 1;
                end
            end
            AnalysisResults.tbExclude(tInd) = CW.CWBias(bInd);
        end

        if RunParams.makePhenoProfileMovie
            drawphenoprofilemovieframe(tInd, WaveResults, ...
                SimResults, SimParams, CW, MFTmodel, DrawParams, Names)
        end

%             if isempty(exclude)
%                 exclude = CW.ncwb - 50;
%             end
%             plot_inds = (exclude(1)-50):10:(exclude(1)+50);
%             plot_inds = plot_inds(plot_inds > 0 & plot_inds < CW.ncwb);
%             subplot(1,2,1)
%             plot(CW.CWBias, squeeze(peak_pheno_pos(:,t)), 'k.', 'LineWidth', lw)
%             xlabel('TB')
%             ylabel('xmax')
%             hold all
%             for p_ind = 1:length(plot_inds)
%                 plot(CW.CWBias(plot_inds(p_ind)), peak_pheno_pos(plot_inds(p_ind),t), 'o', 'LineWidth', lw)
%             end
%             hold off
%             subplot(1,2,2)
%             plot(1000*x, squeeze(SimResults.dynProfile.cellDensity(plot_inds,:,t))', 'LineWidth', lw)
%             xlabel('x')
%             ylabel('cell density')
%             title(['t = ', num2str(t*outDt/60), ' min'])
%             pause
    end

    % find escape flux for each phenotype
    for cwInd = cwIndRange
        nPhenotype = AnalysisResults.waveCellsGMax(tInd) * AnalysisResults.PWaveGMax(cwInd,tInd) * CW.dCWBias; % number of cells of this phenotype from maximal gradient
        
        % assuming group speed is the same as peak speed and is the same as chi f' at boundary, find leak rate
        if CW.CWBias(cwInd) > WaveResults.tbMaxAtArgmaxFz
            AnalysisResults.gsFluxBackCells(cwInd,tInd,1) = (WaveResults.chiMinAtArgmaxFz - MFTmodel.chi(cwInd))/WaveResults.chiMinAtArgmaxFz^2 * WaveResults.c^2 * nPhenotype;
        end
        if CW.CWBias(cwInd) > WaveResults.tbMaxAtChiMin
            AnalysisResults.gsFluxBackCells(cwInd,tInd,2) = (WaveResults.chiMin - MFTmodel.chi(cwInd))/WaveResults.chiMin^2 * WaveResults.c^2 * nPhenotype;
        end
        if CW.CWBias(cwInd) > WaveResults.tbMaxAtAKiA
            AnalysisResults.gsFluxBackCells(cwInd,tInd,3) = (WaveResults.chiMinAtAKiA - MFTmodel.chi(cwInd))/WaveResults.chiMinAtAKiA^2 * WaveResults.c^2 * nPhenotype;
        end
        if CW.CWBias(cwInd) > WaveResults.tbMaxAtAKiAKA
            AnalysisResults.gsFluxBackCells(cwInd,tInd,4) = (WaveResults.chiMinAtAKiAKA - MFTmodel.chi(cwInd))/WaveResults.chiMinAtAKiAKA^2 * WaveResults.c^2 * nPhenotype;
        end
        if CW.CWBias(cwInd) > WaveResults.tbMaxAtF2
            AnalysisResults.gsFluxBackCells(cwInd,tInd,5) = (WaveResults.chiMinAtF2 - MFTmodel.chi(cwInd))/WaveResults.chiMinAtF2^2 * WaveResults.c^2 * nPhenotype;
        end
        if CW.CWBias(cwInd) > WaveResults.tbMaxAtLastEps
            AnalysisResults.gsFluxBackCells(cwInd,tInd,6) = (WaveResults.chiMinAtLastEps - MFTmodel.chi(cwInd))/WaveResults.chiMinAtLastEps^2 * (AnalysisResults.waveCellsGMax(tInd)/SimParams.OD2CellPerdx*SimParams.aA/SimParams.asp)^2 * nPhenotype;
        end
        if CW.CWBias(cwInd) > tbMaxAtChiMinPred
            AnalysisResults.gsFluxBackCells(cwInd,tInd,7) = (chiMinPred - MFTmodel.chi(cwInd))/chiMinPred^2 * (WavePrediction.cPred)^2 * nPhenotype;
        end
        
        % find potential
        chif = MFTmodel.chi(cwInd)*WaveResults.gradient(WaveResults.phenotypeMaxInd(cwInd)); % chi * f' at peak density of this phenotype
        AnalysisResults.cfDriftPhePeak(cwInd,tInd) = chif;
%         AnalysisResults.cpFluxBackCells(cwInd,tInd) = WaveResults.c * ...
%             WaveResults.dynDen(cwInd,1) * SimParams.OD2Cell / SimParams.dx; % um/s * OD * cell/OD / um, number flux
%         AnalysisResults.fpFluxBackCells(cwInd,tInd) = chif * ...
%             WaveResults.dynDen(cwInd,1) * SimParams.OD2Cell / SimParams.dx;
        if isnan(WaveResults.c)
            continue
        end
%         Ui = WaveResults.c*WaveResults.z - MFTmodel.chi(cwInd)*WaveResults.f;
%         Ui = chif*WaveResults.z - MFTmodel.chi(cwInd)*WaveResults.f;
        Ui = AnalysisResults.dxdtMaxGrad(tInd)*WaveResults.z - MFTmodel.chi(cwInd)*WaveResults.f;
        AnalysisResults.U{tInd}(cwInd,:) = Ui;
        
        % draw zeroth order solution
        if RunParams.reDrawZerothOrderSol
            drawzerothordersolution(tInd, cwInd, Ui, WaveResults, SimParams, CW, MFTmodel, DrawParams, Names)
        end
        
        % calculate rate of leaving potential well
        criticalPoints = find(diff(sign(diff(Ui))) ~= 0);
        if numel(criticalPoints) < 2 % no energy well and barrier
            continue
        else
            mu = MFTmodel.mu(cwInd);
            indB = criticalPoints(1) + 1; % index along WaveResults.z at zb
            indM = criticalPoints(end) + 1; % index along WaveResults.z at zm
            d2Uidx2 = [0 diff(Ui, 2) 0]/(SimParams.dx^2);
            d3Uidx3 = ([0 diff(Ui, 3) 0 0]+[0 0 diff(Ui, 3) 0])/2/(SimParams.dx^3);
            d4Uidx4 = [0 0 diff(Ui, 4) 0 0]/(SimParams.dx^4);
            er = sqrt( abs(d2Uidx2(indB)) * d2Uidx2(indM) ) / (2*pi*exp( (Ui(indB)-Ui(indM))/mu )); % escape rate (Risken 5.111)
            imEr = er*(1-mu*(1/8*(d4Uidx4(indB)/d2Uidx2(indB)^2-d4Uidx4(indM)/d2Uidx2(indM)^2)...
                +5/24*(d3Uidx3(indB)^2/d2Uidx2(indB)^3-d3Uidx3(indM)^2/d2Uidx2(indM)^3)...
                )); % improved escape rate (Risken 5.112)
            acEr = mu/(sum(exp(-Ui(indB:end)/mu))*SimParams.dx * ...
                sum(exp(Ui(1:indM)/mu))*SimParams.dx); % more accurate escape rate (Risken (5.109))
            AnalysisResults.fluxBackCells(cwInd,tInd) = er * nPhenotype; % cell/s, number flux
            AnalysisResults.imFluxBackCells(cwInd,tInd) = imEr * nPhenotype;
            AnalysisResults.acFluxBackCells(cwInd,tInd) = acEr * nPhenotype;
            
%             dsd = diff(sign( diff(WaveResults.dynDen(cwInd,:)) ));
%             denMinInd = find(dsd == 2) + 1; % find where first derivative = 0 and second derivative > 0
%             if isempty(denMinInd)
%                 denMinInd = 1;
%             end
%             denMinInd = denMinInd(1); % the first local minimum is real
%             AnalysisResults.cqFluxBackCells(cwInd,tInd) = WaveResults.c * ...
%                 WaveResults.phenotypeMin(cwInd) * SimParams.OD2Cell / SimParams.dx;
%             AnalysisResults.fqFluxBackCells(cwInd,tInd) = chif * ...
%                 WaveResults.phenotypeMin(cwInd) * SimParams.OD2Cell / SimParams.dx;
        end
    end
end
AnalysisResults.peakPhenoSpd(:, AnalysisParams.calculatePeakSpacing+1:end) = (AnalysisResults.peakPhenoPos(:, AnalysisParams.calculatePeakSpacing+1:end) - AnalysisResults.peakPhenoPos(:, 1:end-AnalysisParams.calculatePeakSpacing))/AnalysisParams.calculatePeakSpacing/SimParams.outDt;
AnalysisResults.relSpeedMaxGrad = AnalysisResults.cfDriftMaxGrad + AnalysisResults.mrDriftMaxGrad - ones(CW.ncwb,SimParams.nT) .* AnalysisResults.dxdtMaxGrad;
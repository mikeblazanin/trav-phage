function NumberDensities = leakagedynamics(tInd, AnalysisResults, SimParams, CW, MFTmodel)
% This function generates predicted and simulated number densities over
% time, using the wave profile at some initial time. The algorithm is
% derived in leakage_dynamics.nb.

% OUTPUT_FIELDS
% tInd           : time index to initialize number density prediction
% dCWBias        : dTB used in integration
% simulatedND    : CW.ncwb * SimParams.nT array of simulated number densities over time
% predictedND    : CW.ncwb * SimParams.nT array of predicted number densities over time

% initialization
disp('Calculating leakage dynamics ...')
NumberDensities = [];
NumberDensities.tInd = tInd;
NumberDensities.dCWBias = CW.dCWBias;
NumberDensities.simulatedND = AnalysisResults.PWaveAll.*AnalysisResults.waveCellsAll;
NumberDensities.predictedND = nan(CW.ncwb,SimParams.nT);
NumberDensities.predictedND(:, tInd) = NumberDensities.simulatedND(:, tInd);
indInWave = find(NumberDensities.predictedND(:, tInd) > 0);
bdInd = indInWave(end);
ldFactor = SimParams.KiA * SimParams.aA^2 / SimParams.asp^3; % constant multiplicative factor in leakage dynamics

for t_ind = (tInd+1):SimParams.nT
    NumberDensities.predictedND(:, t_ind) = NumberDensities.predictedND(:, t_ind-1);
    cumDt = 0; % cumulative time passed for the discrete phenotype by phenotype leaking calculation
    while cumDt < SimParams.outDt % leak away the boundary phenotype
        DeltaN = NumberDensities.predictedND(bdInd, t_ind-1) * NumberDensities.dCWBias;
        bdChi = MFTmodel.chi(bdInd);
        bdMu = MFTmodel.mu(bdInd);
        Nother = sum(NumberDensities.predictedND(1:bdInd-1, t_ind-1)) * NumberDensities.dCWBias;
        DeltaT = DeltaN / (ldFactor * bdChi/(bdChi-bdMu)^2 * Nother^3);
        if cumDt + DeltaT <= SimParams.outDt
            NumberDensities.predictedND(bdInd, t_ind) = 0;
            bdInd = bdInd - 1;
        else % leak only a fraction of the last phenotype
            leakFraction = (SimParams.outDt - cumDt) / DeltaT; % guaranteed to be in the open interval (0,1)
            NumberDensities.predictedND(bdInd, t_ind) = (1 - leakFraction) * NumberDensities.predictedND(bdInd, t_ind);
        end
        cumDt = cumDt + DeltaT;
    end
end
end
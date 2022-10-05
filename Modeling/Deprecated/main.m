% error
for asp = 100 % uM , total aspartate concentration
for oxy = 0 % uM, external oxygen level
for ser = 0 % uM, total serine concentration
for nut = 0 % mg/mL, total nutrient concentration
for phi = 0 % proportion of the second population
    %% initialization
    % run script choices
    RunParams = [];
    RunParams.reSimulate            = 0;
    RunParams.reDrawEvolution       = 0;
    RunParams.reDrawSample          = 0;
    RunParams.reAnalyze             = 0;
    RunParams.reDrawZerothOrderSol  = 0;
    RunParams.makePhenoProfileMovie = 0;
    RunParams.reDrawPeakProfile     = 0;
    RunParams.reDrawTBDistribution  = 0;
    RunParams.reDrawLeakage         = 1;
    RunParams.reDrawGradMaxDrift    = 0;
    RunParams.reDrawLeakageDynamics = 0;

    SimParams = [];
    SimParams.concType = 4;
    SimParams.ptbType = 6;
    SimParams.asp = asp;
    SimParams.oxy = oxy;
    SimParams.ser = ser;
    SimParams.nut = nut;
    SimParams.glu = 0;
    SimParams.secondPop = []; % parameters of the second population
    SimParams.secondPop.phi = phi;
    
    [SimParams, DrawParams, CW, MFTmodel, Names] = initializesimulation(SimParams);
    
    if RunParams.reSimulate
        if ~exist(Names.dataDir, 'dir')
            mkdir(Names.dataDir)
        end
        if ~exist(Names.analysisDir, 'dir')
            mkdir(Names.analysisDir)
        end
        %% run simulation
        tic;
        disp(['Running simulation ', Names.simName, '...'])
        SimResults = simulatewave(SimParams, CW, MFTmodel, Names, DrawParams);
        toc
        %% save
        disp(['Saving data file ', Names.simName, '...'])
        save(Names.simName, 'SimParams', 'CW', 'MFTmodel', 'Names', 'SimResults', '-v7.3')
    else
        if ~exist(Names.simName, 'file')
            disp(['The simulation ', Names.simName, ' does not exist, skipping...'])
            continue
        end
        
        % load
        disp(['Loading data file ', Names.simName, '...'])
        load(Names.simName)
    end
    
    %% initialize plot variables
    [AnalysisParams, DrawParams] = initializeanalysis(SimParams, SimResults, DrawParams, Names);
    
    %% draw wave profile evolution
    if RunParams.reDrawEvolution
        drawprofileevolution(AnalysisParams, SimResults, SimParams, DrawParams, Names)
    end
    
    %% loop over a few times
    if RunParams.reDrawSample
        disp('Drawing profiles at sampled times...')
        for plotTimeInd = 1:DrawParams.nTInds
            tInd = DrawParams.tInds(plotTimeInd);
            disp(['t = ', num2str(tInd * SimParams.outDt), 's'])
            WaveResults = analyzewavesnapshot(tInd, AnalysisParams.wb(plotTimeInd), AnalysisParams.wf(plotTimeInd), ...
                        AnalysisParams, SimResults, SimParams, CW, MFTmodel);
            if isempty(WaveResults.dynDen)
                disp(['Wave does not exist at t = ', num2str(tInd * SimParams.outDt), ' s, skipping...'])
                continue
            end
            drawsample(tInd, WaveResults, SimParams, CW, MFTmodel, DrawParams, Names)
            drawshapeprediction(tInd, WaveResults, SimParams, CW, MFTmodel, DrawParams, Names)
        end
    end
    
    % calcualte wave profile time series
    if RunParams.reAnalyze
        disp('Analyzing wave peaks...')
        AnalysisResults = analyzewavedynamics(AnalysisParams, SimResults, ...
            SimParams, CW, MFTmodel, DrawParams, Names, RunParams);
        %% save
        if ~RunParams.reDrawZerothOrderSol % don't save partial results produced when this option is on
            disp(['Saving analysis file ', Names.analysisName, '...'])
            save(Names.analysisName, 'AnalysisResults', 'AnalysisParams', 'DrawParams', '-v7.3')
        else
            if exist(Names.analysisName, 'file') % get back the full AnalysisResults
                disp(['Loading analysis file ', Names.analysisName, '...'])
                load(Names.analysisName)
            else % if no file, clean up partial results
                clear AnalysisResults
            end
        end
    end
    
    %% make sure AnalysisResults exists
    if exist(Names.analysisName, 'file')
        disp(['Loading analysis file ', Names.analysisName, '...'])
        load(Names.analysisName)
    else
        disp(['The analysis file ', Names.analysisName, ' does not exist, skipping ...'])
        continue
    end
    
    %% make movie
    if RunParams.makePhenoProfileMovie
        writephenoprofilemovie(SimParams, Names)
    end
    
    %% plot multiple
    if ~exist('AnalysisResultsMultiple', 'var')
        AnalysisResultsMultiple = [];
    end
    nCases = 1;
    nCase = 1;
    AnalysisResultsMultiple.peakSpdMultiple{nCase}   = AnalysisResults.peakSpd;
    AnalysisResultsMultiple.peakDenMultiple{nCase}   = AnalysisResults.peakDen;
    AnalysisResultsMultiple.waveCellsMultiple{nCase} = AnalysisResults.waveCellsAll;
    AnalysisResultsMultiple.tbExcludeMultiple{nCase} = AnalysisResults.tbExclude;
    
    %% show profiles
    if RunParams.reDrawPeakProfile
        plotMultiple = false;
        drawpeakprofile(AnalysisResults, SimParams, DrawParams, Names, ...
            plotMultiple, nCases, AnalysisResultsMultiple)
    end
    
    %% show tumble bias distribution
    if RunParams.reDrawTBDistribution
        drawtbdistribution(AnalysisResults, SimParams, CW, DrawParams, Names)
    end

    %% show cell flux at the back
    if RunParams.reDrawLeakage
        plotMultiple = false;
        drawleakage(AnalysisParams, AnalysisResults, SimParams, CW, MFTmodel, DrawParams, Names, plotMultiple)
    end
    
    %% show cell flux at gradient max
    if RunParams.reDrawGradMaxDrift
        drawgradmaxdrift(AnalysisResults, SimParams, CW, DrawParams, Names)
    end
    
    %% show leakage dynamics
    if RunParams.reDrawLeakageDynamics
        NumberDensities = leakagedynamics(80, AnalysisResults, SimParams, CW, MFTmodel);
        drawleakagedynamics(NumberDensities, SimParams, CW, DrawParams, Names)
    end
end
end
end
end
end
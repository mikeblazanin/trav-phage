function [AnalysisParams, DrawParams] = initializeanalysis(SimParams, SimResults, DrawParams, Names)
% This file initialize analysis parameters.
AnalysisParams = [];

DrawParams.tInds = 4:4:29; % time index at which to plot sample outputs
DrawParams.nTInds = numel(DrawParams.tInds);
DrawParams.waveBound = [-3000 1500]; % um, boundary of plots in moving coordinate
DrawParams.cMax = 20; % ylim when plotting wave speed

AnalysisParams.calculatePeakSpacing = 3; % time index spacing to calculate peak speed

% The analysis only works for a single population and no attempt has been
% made to adopt to two populations cases.
if ndims(SimResults.dynProfile.cellDensity) > 3
    disp('The rest of the code does not deal with two populations.')
    return
end

if exist(Names.analysisName, 'file')
    disp(['Loading analysis file ', Names.analysisName, ' ...'])
    load(Names.analysisName)
else % determine where the wave is, this part relies on handpicking
    disp(['The analysis file ', Names.analysisName, ' does not exist, requiring handpicking...'])
    figure
    plot(squeeze(sum(SimResults.dynProfile.cellDensity(:,:,DrawParams.tInds))))
    AnalysisParams.wo = -20+25*(1:DrawParams.nTInds); % array of absolute index along x where the wave starts
    AnalysisParams.wb = 50+5*(1:DrawParams.nTInds); % array of index length from the back to the peak of the wave
    AnalysisParams.wf = 20+5*(1:DrawParams.nTInds); % array of index length fomr the front to the peak of the wave
    
    % use peaks
    if SimParams.concType == 4
        AnalysisParams.wo = [];
        AnalysisParams.wb = [];
        AnalysisParams.wf = [];
        
        for i = 1:length(DrawParams.tInds)
            [~,locso] = findpeaks(squeeze(sum(SimResults.dynProfile.cellDensity(:,:,DrawParams.tInds(i)),1)));
            [~,locsb] = findpeaks(-squeeze(sum(SimResults.dynProfile.cellDensity(:,:,DrawParams.tInds(i)),1)));
            if ~isempty(locso)
                locsb = locsb(locsb<locso(end));
            end
            if isempty(locsb)
                locsb = 1;
            end
            if isempty(locso)
                locso = 1;
            end

            AnalysisParams.wo(DrawParams.tInds(i)) = locso(end);
            AnalysisParams.wb(DrawParams.tInds(i)) = locso(end)-locsb(end);
%             AnalysisParams.wf(DrawParams.tInds(i)) = size(SimResults.dynProfile.cellDensity,2) - locso(end);
            maxDens = max(squeeze(sum(SimResults.dynProfile.cellDensity(:,:,DrawParams.tInds(i)),1))); 
            AnalysisParams.wf(DrawParams.tInds(i)) = find(squeeze(sum(SimResults.dynProfile.cellDensity(:,locso(end):end,DrawParams.tInds(i)),1))<1e-6*maxDens,1,'first');
            
            plot(locso(end),squeeze(sum(SimResults.dynProfile.cellDensity(:,locso(end),DrawParams.tInds(i)),1)),'ko')
            plot(locsb(end),squeeze(sum(SimResults.dynProfile.cellDensity(:,locsb(end),DrawParams.tInds(i)),1)),'ro')
        end
       
    end
end
end
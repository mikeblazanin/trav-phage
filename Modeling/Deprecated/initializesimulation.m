function [SimParams, DrawParams, CW, MFTmodel, Names] = initializesimulation(SimParams)

disp('Initializing simulation ...')
% figure parameterss
DrawParams = [];
DrawParams.rez  = 300; % save figure resolution
DrawParams.figW = 5; % inch, figure width
DrawParams.lw   = 3; % line width

Names = [];

% second population
if ~isfield(SimParams, 'secondPop')
    SimParams.secondPop = []; % parameters of the second population
end
if ~isfield(SimParams.secondPop, 'phi')
    SimParams.secondPop.phi = 0; % proportion of the second population
end
if SimParams.secondPop.phi > 0
    SimParams.secondPop.speedRatio = 1; % ratio of the speed of the second population
    SimParams.secondPop.aARatio = 0.5; % ratio of the aspartate consumption of the second population
    SimParams.secondPop.aORatio = 1; % ratio of the oxygen consumption of the second population
    SimParams.secondPop.aSRatio = 2; % ratio of the serine consumption of the second population
    SimParams.secondPop.aNRatio = 1; % ratio of the nutrient consumption of the second population
    SimParams.secondPop.sGRatio = 1; % ratio of the glutamate secretion of the second population
    SimParams.secondPop.cGRatio = 1; % ratio of the glutamate consumption of the second population
end

% simulate asp
if ~isfield(SimParams, 'asp')
    SimParams.asp = 200; % uM , total aspartate concentration
end
if SimParams.asp > 0
    Names.aStr = 'Asp_';
else
    Names.aStr = '';
end

% simulate oxy
if ~isfield(SimParams, 'oxy')
    SimParams.oxy = 250; % uM, external oxygen level
end
if SimParams.oxy > 0
    Names.oStr = 'Oxy_';
else
    Names.oStr = '';
end

% simulate ser
if ~isfield(SimParams, 'ser')
    SimParams.ser = 0; % uM, total serine concentration
end
if SimParams.ser > 0
    Names.sStr = 'Ser_';
else
    Names.sStr = '';
end

% simulate nutrient (non-attractant)
if ~isfield(SimParams, 'nut')
    SimParams.nut = 0.5; % mg/mL, total nutrient concentration
end
if SimParams.nut > 0
    Names.nStr = 'Nut_';
else
    Names.nStr = '';
end

% simulate glu
if ~isfield(SimParams, 'glu')
    SimParams.glu = 0; % whether to simulate glutamate as a self-attractant
end
if SimParams.glu > 0
    Names.gStr = 'Glu_';
else
    Names.gStr = '';
end

% initial concentration profile of asp, oxy, ser, and nutrient
if ~isfield(SimParams, 'concType')
    SimParams.concType = 0; % initial concentration profile type
                            % 0 = Hill-shaped
                            % 1 = flat
                            % 2 = centrifuge simulated
                            % 3 = theoretical prediction of traveling wave
                            % 4 = initialize wave using time/length scale analyses
end
switch SimParams.concType
    case 0
        Names.cStr = 'Hill';
    case 1
        Names.cStr = 'Flat';
    case 2
        Names.cStr = 'Cent';
    case 3
        Names.cStr = 'Theo';
    case 4
        Names.cStr = 'Scal';
end
Names.cStr = [Names.cStr, 'Conc_'];

% initial tumble bias distribution
if ~isfield(SimParams, 'ptbType')
    SimParams.ptbType = 4; % initial tumble bias distribution
                           % 0 = Bell-shaped
                           % 1 = uniform
                           % 2 = custom distribution
                           % 3 = Adam's WT distribution
                           % 4 = Xiongfei's interpolation
                           % 5 = Sum of 2 hat distributions
                           % 6 = log-normal
                           % 7 = Gamma distribution
end
SimParams.ptbParams = [];
switch SimParams.ptbType
    case 0
        Names.pStr = 'Bell';
    case 1
        Names.pStr = 'Unif';
    case 2
        Names.pStr = 'Cstm';
    case 3
        Names.pStr = 'Adam';
    case 4
        Names.pStr = 'XFei';
    case 5
        Names.pStr = 'Hat2';
        SimParams.ptbParams.hatWeight1 = 1;
        SimParams.ptbParams.hatLeft1   = 0;
        SimParams.ptbParams.hatRight1  = 0.6;
        SimParams.ptbParams.hatWeight2 = 0;
        SimParams.ptbParams.hatLeft2   = 0;
        SimParams.ptbParams.hatRight2  = 1;
    case 6
        Names.pStr = 'LogN'; % default mean = 0.3
        % mean 0.3 std 0.01
%         SimParams.ptbParams.sigma = sqrt(log((0.01/0.3)^2 + 1));
%         SimParams.ptbParams.mu = log(0.3) - SimParams.ptbParams.sigma^2/2;
        % mean 0.3 std 0.1
        SimParams.ptbParams.sigma = 0.3246;
        SimParams.ptbParams.mu = -1.2567;
        % mean 0.3 std 0.26
%         SimParams.ptbParams.sigma = 0.7485;
%         SimParams.ptbParams.mu = -1.4841;
        % mean 0.5 std 0.01
%         SimParams.ptbParams.sigma = 0.02;
%         SimParams.ptbParams.mu = -0.6933;
        % mean 0.5 std 0.1
%         SimParams.ptbParams.sigma = 0.1980;
%         SimParams.ptbParams.mu = -0.7128;
        % mean 0.6 std 0.1
%         SimParams.ptbParams.sigma = 0.1655;
%         SimParams.ptbParams.mu = -0.5245;
    case 7
        Names.pStr = 'Gama';
%             SimParams.ptbParams.a = 30;
%             SimParams.ptbParams.b = 0.3 / SimParams.ptbParams.a;
%             SimParams.ptbParams.a = 3;
%             SimParams.ptbParams.b = 0.10087;
        SimParams.ptbParams.a = 1;
        SimParams.ptbParams.b = 0.37637;
end
Names.pStr = [Names.pStr, 'TB_'];

% receptor nonlinearity
if ~isfield(SimParams, 'recType')
    SimParams.recType = 0; % receptor nonlinearity
                           % 0 = nonlinear
                           % 1 = log
                           % 2 = linear
end
switch SimParams.recType
    case 0
        Names.rStr = 'NL';
    case 1
        Names.rStr = 'LG';
    case 2
        Names.rStr = 'LI';
end
Names.rStr = [Names.rStr, 'rec_'];

% dimensionality
SimParams.dim = 1;
switch SimParams.dim
    case 1
        Names.dStr = '1D';
    case 2
        Names.dStr = '2D';
end

% aspartate parameters
SimParams.DA  = 500; % um^2/s, aspartate diffusion coefficient, A method for the determination of diffusion coefficients for small molecules in aqueous solution, Analytical Biochemistry, Volume 166, Issue 2, 1 November 1987, Pages 335-341
SimParams.aA  = 9*0.84/60; % uM/OD/s, asp consumption rate when cell density is low (oxygen is maximal), 1e-12 umol/cell/min = 0.84/60 uM/OD/s
SimParams.KA  = 0.5; % uM, Michaelis-Menten constant in asp consumption
SimParams.lAO = 0.27; % fraction of asp consumption without oxygen
SimParams.KAO = 250; % uM, scaling constant in oxygen-dependence of asp consumption
SimParams.HAO = 1.1; % steepness in oxygen-dependence of asp consumption
SimParams.KiA  = 3.5;  % uM, association constant of asp with chemoreceptor
SimParams.KaA  = 55000; % uM, dissociation constant of asp with chemoreceptor

% normalize consumption when oxygen is not present
if SimParams.oxy == 0
    lowDensityWeight = 0.; % hyperparameter used to ensure comparison of cases with and without oxygen, try 0.9
    SimParams.aA   = SimParams.aA * (SimParams.lAO*lowDensityWeight + 1 - lowDensityWeight);
    SimParams.lAO  = 1;
end

% oxygen parameters
if SimParams.oxy > 0
    SimParams.DO = 2500; % um^2/s, oxygen diffusion coefficient, http://compost.css.cornell.edu/oxygen/oxygen.diff.water.html
    SimParams.aO = 70*0.84/60; % uM/OD/s, maximal oxygen consumption rate, 1e-12 umol/cell/min = 0.84/60 uM/OD/s
    SimParams.KO = 1; % uM, Michaelis-Menten constant in oxygen consumption
    SimParams.kO = 0.02; % 1/s, oxygen transfer rate from PDMS to liquid interface
    SimParams.KiO = SimParams.KiA; % uM, association constant of oxy with chemoreceptor
    SimParams.KaO = SimParams.KaA; % uM, dissociation constant of oxy with chemoreceptor
    SimParams.chiO = 0; % chi_O / chi_A
end

% serine parameters
if SimParams.ser > 0
    SimParams.DS = 500; % um^2/s, serine diffusion coefficient
    SimParams.aS = 7.5/60; % uM/OD/s, maximal serine consumption rate
    SimParams.KS = 0.5; % uM, Michaelis-Menten constant in serine consumption
    SimParams.KiS = 0.9/0.3*SimParams.KiA; % uM, association constant of glu with chemoreceptor
    SimParams.KaS = 0.9/0.3*SimParams.KaA; % uM, dissociation constant of glu with chemoreceptor
    SimParams.chiS = 1; % chi_S / chi_A
end

% nutrient parameters
if SimParams.nut > 0
    SimParams.DN = 500; % um^2/s, nutrient diffusion coefficient
    SimParams.aN = 0.05; % mg/mL/OD/s, nutrient consumption rate
    SimParams.KN = 0.001; % mg/mL, Michaelis-Menten constant in nutrient consumption
end

% glutamate parameters
if SimParams.glu > 0
    SimParams.DG = 500; % um^2/s, glutamate diffusion coefficient
    SimParams.sG = 1; % uM/OD/s, glutamate secretion rate, 1 uM/OD/s = 7.17e5 molecules/cells/s
    SimParams.dG = 5e-3; % 1/s, glutamate degradataion rate
    SimParams.cG = 1.4e-5; % 1/OD/s, glutamate consumption rate
    SimParams.KiG = 50/0.3*SimParams.KiA; % uM, association constant of glu with chemoreceptor
    SimParams.KaG = 50/0.3*SimParams.KaA; % uM, dissociation constant of glu with chemoreceptor
    SimParams.chiG = 1; % chi_G / chi_A
end

% cell model
CW = initializeCW(SimParams);
MFTmodel = initializeMFTmodel(CW, SimParams);

% cell parameters
SimParams.NCells = 5e4; % total number of cells initialized
SimParams.growthRate = 0; % cell growth rate
SimParams.inherit = 0; % how phenotype is being passed on, only for growthRate > 0
SimParams.PNew = repmat(CW.P*CW.dCWBias,1,length(CW.P)); % probability of phenotypes in a new generation of cells

% discretization and channel set up
SimParams.totTime = 30*60; % s, total simulation time
SimParams.outDt   = 10; % s, the output time interval
if ~isfield(SimParams, 'dx')
    SimParams.dx = 20; % um, spatial discretization
end
SimParams.dt      = SimParams.dx^2/4000; % s, relation with dx to ensure stability
SimParams.channelLength = 30000; % um, relationship with asp is empirical
SimParams.OD1           = 8.4e8; % cells/ml in OD1 
SimParams.width         = 600; % um, channel width
SimParams.height        = 14; % um, channel height
SimParams.OD2CellPerdx  = SimParams.width*SimParams.height*SimParams.OD1/1e12; % area * density = number of cells per OD per length in um
SimParams.V             = 1740; % um, initial profile length scale
SimParams.H             = 3; % initial profile steepness
% scale the mesh and chooses nice initial condition
if SimParams.concType == 4
    % recompute grid
    c = SimParams.NCells/SimParams.OD2CellPerdx*SimParams.aA/SimParams.asp;
    meanChi = sum(CW.P*CW.dCWBias.*MFTmodel.chi); % um^2/s
%     stdChi = sqrt(sum(CW.P*CW.dCWBias.*(MFTmodel.chi-meanChi).^2)); % um^2/s
    SimParams.tScale = meanChi/c^2; % s, time for wave to move one characteristic length
    SimParams.xScale = meanChi/c; % um, characteristic length of wave
    SimParams.rScale = SimParams.NCells/SimParams.OD2CellPerdx/SimParams.xScale; % OD
    
    SimParams.scaleSims = true;
    if SimParams.scaleSims
        SimParams.dx = round(SimParams.xScale/10,3);
        SimParams.dt = round(SimParams.dx^2/meanChi/10,3);
        SimParams.outDt = round(SimParams.tScale); % max(SimParams.outDt,SimParams.dt);
        SimParams.totTime = 30*SimParams.tScale;
        SimParams.channelLength = 20*SimParams.xScale;
    end
end
% continue set up
SimParams.OD2Cell    = SimParams.OD2CellPerdx * SimParams.dx; % number of cells per OD per simulation block
SimParams.nT         = floor(SimParams.totTime/SimParams.outDt); % total number of output times
SimParams.time       = 0:SimParams.outDt:(SimParams.nT-1)*SimParams.outDt; % array of output times
SimParams.L          = floor(SimParams.channelLength/SimParams.dx); % number of spatial index
SimParams.Drho       = MFTmodel.mu*ones(1,SimParams.L); % cell spatial diffusion coefficient
SimParams.x          = 0:SimParams.dx:((SimParams.L-1)*SimParams.dx); % um, position along channel
SimParams.centrifuge = 0; % whether to include centrifugal force
if SimParams.centrifuge == 1 % centrifuge parameters
    SimParams.r0 = 70000; % distance in um from the beginning of the channel to the center of the SimParams.centrifuge, 7 cm
    SimParams.w0 = 3000*2*pi/60; % spinning speed, 3000 rpm
    SimParams.sedCoeff = 1e3*(1e-6)^3/(6*pi*1e-3*1e-6); % sedimentation coefficient of E. coli, m/(6*pi*eta*r)
    SimParams.DL = round(10000/SimParams.dx); % cells loaded in the first 10000 um of the channel
else
    SimParams.DL = round(800/SimParams.dx); % cells loaded in the first 800 um of the channel and diffused
end
if SimParams.dim == 2
    SimParams.r = SimParams.x - 1/2*SimParams.dx; % um, radial position in 2D
end

% names
Names.initCondDir = [Names.aStr, Names.oStr, Names.sStr, Names.nStr, Names.gStr, Names.cStr, Names.pStr, Names.rStr, Names.dStr, '\'];
Names.paramDir = [num2str(SimParams.asp), '_', num2str(SimParams.oxy), '_', ...
    num2str(SimParams.ser), '_', num2str(SimParams.nut), '_', ...
    num2str(SimParams.glu), '_', num2str(SimParams.totTime), '_', ...
    num2str(SimParams.outDt), '_', num2str(SimParams.channelLength), '_', ...
    num2str(SimParams.dt), '_', num2str(SimParams.dx), '\'];
Names.subCaseDir = '';
switch SimParams.ptbType
    case 5
        Names.subCaseDir = ['ptb5_', num2str(SimParams.ptbParams.hatWeight1), '_', ...
            num2str(SimParams.ptbParams.hatLeft1), '_', ...
            num2str(SimParams.ptbParams.hatRight1), '_', ...
            num2str(SimParams.ptbParams.hatWeight2), '_', ...
            num2str(SimParams.ptbParams.hatLeft2), '_', ...
            num2str(SimParams.ptbParams.hatRight2), '_', Names.subCaseDir];
    case 6
        Names.subCaseDir = ['ptb6_', num2str(SimParams.ptbParams.mu), '_', num2str(SimParams.ptbParams.sigma), '_', Names.subCaseDir];
    case 7
        Names.subCaseDir = ['ptb7_', num2str(SimParams.ptbParams.a), '_', num2str(SimParams.ptbParams.b), '_', Names.subCaseDir];
end
if SimParams.glu > 0
    Names.subCaseDir = ['sG=', num2str(SimParams.sG), ...
        '_dG=', num2str(SimParams.dG), ...
        '_cG=', num2str(SimParams.cG), '_', Names.subCaseDir];
end
if SimParams.secondPop.phi > 0
    Names.subCaseDir = ['phi=', num2str(SimParams.secondPop.phi), ...
        '_v=', num2str(SimParams.secondPop.speedRatio), ...
        '_aA=', num2str(SimParams.secondPop.aARatio), ...
        '_aO=', num2str(SimParams.secondPop.aORatio), ...
        '_aS=', num2str(SimParams.secondPop.aSRatio), ...
        '_aN=', num2str(SimParams.secondPop.aNRatio), ...
        '_sG=', num2str(SimParams.secondPop.sGRatio), ...
        '_cG=', num2str(SimParams.secondPop.cGRatio), ...
        '_', Names.subCaseDir];
end
if SimParams.growthRate > 0
    Names.subCaseDir = ['growthRate=', num2str(SimParams.growthRate), '_', ...
        Names.subCaseDir];
end

% normalize cell counts at one point in the simulation to diminish the effect of the separation process
SimParams.cellNormalize = 0;
Names.normStr = [];
if SimParams.cellNormalize
    SimParams.normalizeTime = 900; % s, time at which normalization happens
    SimParams.normalizeNCells = 3e4; % number of cells in the wave to normalize to
    SimParams.normalizePosition = 5000; % um, channel position where cell count starts for the wave
    Names.subCaseDir = ['norm_', num2str(SimParams.normalizeTime), '_', ...
        num2str(SimParams.normalizeNCells), '_', ...
        num2str(SimParams.normalizePosition), '_', Names.subCaseDir];
end

if ~isempty(Names.subCaseDir)
	if Names.subCaseDir(end) == '_'
	    Names.subCaseDir = Names.subCaseDir(1:end-1);
	end
	Names.subCaseDir = [Names.subCaseDir, '\'];
end
Names.caseDir = [Names.paramDir, Names.subCaseDir];
Names.rootDir = '..\..\';
Names.dataDir = [Names.rootDir, 'data\SCRATCH\wave\', Names.initCondDir, Names.caseDir];
Names.analysisDir = [Names.rootDir, 'analysis\wave\', Names.initCondDir, Names.caseDir];
Names.phenoProfileMovieDir = [Names.analysisDir, 'pheno_profile\'];
Names.simName = [Names.dataDir, 'results.mat'];
Names.analysisName = [Names.analysisDir, 'analysis.mat'];
end
function CW = initializeCW(SimParams)
% This function initializes CW bias and distribution.

CW.dCWBias = 0.01;
CW.CWBias  = CW.dCWBias:CW.dCWBias:0.99;
% define probability density
switch SimParams.ptbType
    case 0
        % CW.P = lognpdf(CW.CWBias,log(0.295),0.308);
        CW.P = lognpdf(CW.CWBias,log(0.318),0.298);
    case 1 % uniform TB distribution
        CW.P = unifpdf(CW.CWBias,0,1);
    case 2 % custom/empirical distributions
        disp('Loading custom tumble bias distribution...')
        load('custom_TB_distribution.mat') % obtained from get_TB_distribution.m
        CW.P = f;
    case 3 % wild type distribution in Adam's simulation
        disp('Loading Adam''s WT tumble bias distribution...')
        load('Adam_WT_TB_distribution.mat') % obtained from showTBtauScatter.m in Dropbox (emonetlab)\paper_wave_diversity\others\Adam_sims_phenotypes\
        CW.P = f;
    case 4 % Xiongfei's interpolated TB distribution
        CWBiasExp=[0,0.0100000000000000,0.0200000000000000,0.0300000000000000,0.0400000000000000,0.0500000000000000,0.0600000000000000,0.0700000000000000,0.0800000000000000,0.0900000000000000,0.100000000000000,0.110000000000000,0.120000000000000,0.130000000000000,0.140000000000000,0.150000000000000,0.160000000000000,0.170000000000000,0.180000000000000,0.190000000000000,0.200000000000000,0.210000000000000,0.220000000000000,0.230000000000000,0.240000000000000,0.250000000000000,0.260000000000000,0.270000000000000,0.280000000000000,0.290000000000000,0.300000000000000,0.310000000000000,0.320000000000000,0.330000000000000,0.340000000000000,0.350000000000000,0.360000000000000,0.370000000000000,0.380000000000000,0.390000000000000,0.400000000000000,0.410000000000000,0.420000000000000,0.430000000000000,0.440000000000000,0.450000000000000,0.460000000000000,0.470000000000000,0.480000000000000,0.490000000000000,0.500000000000000,0.510000000000000,0.520000000000000,0.530000000000000,0.540000000000000,0.550000000000000,0.560000000000000,0.570000000000000,0.580000000000000,0.590000000000000,0.600000000000000,0.610000000000000,0.620000000000000,0.630000000000000,0.640000000000000,0.650000000000000,0.660000000000000,0.670000000000000,0.680000000000000,0.690000000000000,0.700000000000000];
        PExp=[0,0.00302823490236943,0.00498640060510877,0.00740654246062159,0.0129576724735464,0.0282364445536801,0.0342746634745129,0.0605832532805051,0.0797351468975508,0.105863629784042,0.144957679208212,0.222281447976504,0.281608928805501,0.373987712933409,0.512065751136908,0.691152325282671,0.919556834551468,1.19485305314468,1.51155357679430,1.88448868628358,2.26163802431774,2.67050119556754,3.04502589236555,3.42156794182481,3.78827243915744,4.11858437055020,4.26770099008759,4.46704480752076,4.50323811436766,4.50544966346908,4.42742856648661,4.31612752213579,4.11316467718262,3.94529362322380,3.65109481682462,3.23848249620281,2.98678093597204,2.69627429599337,2.36626462681502,2.10807472740901,1.95992530735338,1.71658120173055,1.53144380527692,1.43826939988074,1.27455605283792,1.14836389445003,1.01135546697596,0.907225653941910,0.780385038989034,0.685196631521439,0.645362573778435,0.620271119730271,0.579787157033761,0.525204939995039,0.523233085565830,0.466362640857505,0.445131198472857,0.399401650650870,0.399447904195710,0.371516048640824,0.380986255044554,0.393076006538824,0.375160618566594,0.357271202286408,0.341376560968248,0.322917046544370,0.291493170711647,0.311004430859554,0.290374229335908,0.272975724047314,0.162176873861474];
%         CW.dCWBias = 0.005;
        CW.CWBias = CW.dCWBias:CW.dCWBias:CWBiasExp(end);
        CW.P = interp1(CWBiasExp,PExp,CW.CWBias,'spline');
    case 5 % sum of 2 hat TB distributions
        CW.P = SimParams.ptbParams.hatWeight1 * unifpdf(CW.CWBias, SimParams.ptbParams.hatLeft1, SimParams.ptbParams.hatRight1)...
            + SimParams.ptbParams.hatWeight2 * unifpdf(CW.CWBias, SimParams.ptbParams.hatLeft2, SimParams.ptbParams.hatRight2);
    case 6 % general log-normal distribution
        CW.P = lognpdf(CW.CWBias, SimParams.ptbParams.mu, SimParams.ptbParams.sigma);
        CW.P = CW.P / sum(CW.P) / CW.dCWBias;
        m = sum(CW.CWBias .* CW.P * CW.dCWBias);
        s = sqrt(sum((CW.CWBias - m).^2 .* CW.P * CW.dCWBias));
        disp(['Choosing log normal distribution with mean ', num2str(m), ' and std ', num2str(s)])
    case 7 % gamma distribution
        CW.P = gampdf(CW.CWBias, SimParams.ptbParams.a, SimParams.ptbParams.b);
        CW.P = CW.P / sum(CW.P) / CW.dCWBias;
        m = sum(CW.CWBias .* CW.P * CW.dCWBias);
        s = sqrt(sum((CW.CWBias - m).^2 .* CW.P * CW.dCWBias));
        disp(['Choosing gamma distribution with mean ', num2str(m), ' and std ', num2str(s)])
end
CW.P = CW.P / sum(CW.P) / CW.dCWBias;
CW.P = CW.P(:);
CW.CWBias = CW.CWBias(:);
CW.ncwb = numel(CW.CWBias);
end
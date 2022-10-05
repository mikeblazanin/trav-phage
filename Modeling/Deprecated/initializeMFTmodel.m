function MFTmodel = initializeMFTmodel(CW, SimParams)
% This function initializes MFT model parameters.

    % motor parameters
    MFTmodel.g1    = 40;
    MFTmodel.g2    = 40;
    MFTmodel.Kd    = 3.06;
    MFTmodel.w     = 1.3;
    MFTmodel.alpha = 5.1138;
    MFTmodel.N     = 6;
    if SimParams.ptbType == 4
        MFTmodel.g1    = 28;
        MFTmodel.g2    = 28;
        MFTmodel.Kd    = 2.9;
        MFTmodel.w     = 3.8;
    end
    
    % motor functions
    MFTmodel.Yp      = MFTmodel.g2*MFTmodel.Kd./(MFTmodel.g2-MFTmodel.g1/2+log(1./CW.CWBias-1))-MFTmodel.Kd;
    MFTmodel.deltaG  = MFTmodel.g1/4-MFTmodel.g2/2*(MFTmodel.Yp./(MFTmodel.Kd+MFTmodel.Yp));
    MFTmodel.lambdaT = MFTmodel.w*exp(MFTmodel.deltaG);
    MFTmodel.lambdaR = MFTmodel.w*exp(-MFTmodel.deltaG);
    MFTmodel.CWBias  = ((MFTmodel.lambdaR)./((MFTmodel.lambdaT)+(MFTmodel.lambdaR)));
    MFTmodel.a       = MFTmodel.Yp/MFTmodel.alpha;

    % physical parameters
    MFTmodel.Drot  = 0.062*2;
    MFTmodel.v     = 28*4/pi; %2D run speed 26um/s
    MFTmodel.theta = 1-0.1564;

    % MFT results
    MFTmodel.lambdaRdiff = -(MFTmodel.lambdaR).*MFTmodel.a.*(MFTmodel.a-1) ...
        .*(MFTmodel.alpha*MFTmodel.g2*MFTmodel.Kd)./(2*(MFTmodel.Kd+MFTmodel.alpha*MFTmodel.a).^2);
    if SimParams.ptbType == 4
%         cwBin = [0 0.0250000000000000,0.0750000000000000,0.125000000000000,0.175000000000000,0.225000000000000,0.275000000000000,0.325000000000000,0.375000000000000,0.425000000000000,0.475000000000000,0.525000000000000,0.575000000000000,0.625000000000000,0.675000000000000,0.725000000000000];
%         muExp = [1100.37707662951,1060.37707662951,546.954659772813,376.328326298169,268.359924863403,194.401038526061,143.168042271451,107.493605443718,82.6199648547623,63.9732110490859,49.7344586006807,38.8063536234783,31.0854022068343,23.8371937085902,17.5907565702979,19.5435281626708];
%         MFTmodel.mu       = interp1(cwBin, muExp, CW.CWBias, 'spline');
        MFTmodel.mu       = MFTmodel.v^2*(1-MFTmodel.CWBias)./3 ...
            ./(MFTmodel.lambdaR*MFTmodel.theta+MFTmodel.Drot); % spatial diffusion coefficient
%         MFTmodel.tau      = 3;
%         MFTmodel.tau      = 15*exp(-2.3*MFTmodel.CWBias); % Xiongfei's fit
%         MFTmodel.chi      = MFTmodel.mu.*MFTmodel.theta*MFTmodel.N*MFTmodel.tau ...
%                 ./(1+MFTmodel.tau.*(MFTmodel.lambdaR*MFTmodel.theta+MFTmodel.Drot));
%         MFTmodel.chi      = MFTmodel.theta*MFTmodel.v^2*MFTmodel.N*MFTmodel.tau ...
%             .*(MFTmodel.lambdaRdiff./(1+MFTmodel.lambdaR./MFTmodel.lambdaT) ...
%             ./(MFTmodel.lambdaR*MFTmodel.theta+MFTmodel.Drot)) ...
%             ./(1+MFTmodel.tau.*(MFTmodel.lambdaR*MFTmodel.theta+MFTmodel.Drot))/3;
        MFTmodel.chi      = MFTmodel.mu*22./(1+(0.01./MFTmodel.CWBias).^5).*(1+1./(1+MFTmodel.CWBias/0.01))./(1+0.15*MFTmodel.CWBias.^2);
    else
        MFTmodel.mu       = MFTmodel.v^2*(1-MFTmodel.CWBias)./3 ...
            ./(MFTmodel.lambdaR*MFTmodel.theta+MFTmodel.Drot); % spatial diffusion coefficient
        MFTmodel.tau      = 20*exp(-6*MFTmodel.CWBias);
        MFTmodel.chi      = MFTmodel.theta*MFTmodel.v^2*MFTmodel.N*MFTmodel.tau ...
            .*(MFTmodel.lambdaRdiff./(1+MFTmodel.lambdaR./MFTmodel.lambdaT) ...
            ./(MFTmodel.lambdaR*MFTmodel.theta+MFTmodel.Drot)) ...
            ./(1+MFTmodel.tau.*(MFTmodel.lambdaR*MFTmodel.theta+MFTmodel.Drot))/3;
%         MFTmodel.chi      = MFTmodel.mu;
    end
end
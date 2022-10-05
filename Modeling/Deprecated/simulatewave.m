function SimResults = simulatewave(SimParams, CW, MFTmodel, Names, DrawParams)    
% This function calculates the chemotactic coefficients based on the PLOS Comp.

%% initialization
disp('Starting wave simulation ...')
figure('visible','off')
subplot(3,1,1)
plot(MFTmodel.CWBias, MFTmodel.chi)
xlabel('TB')
ylabel('\chi (\mum^2/s)')
subplot(3,1,2)
plot(MFTmodel.CWBias, MFTmodel.mu)
xlabel('TB')
ylabel('\mu (\mum^2/s)')
subplot(3,1,3)
plot(MFTmodel.CWBias, MFTmodel.chi./MFTmodel.mu)
xlabel('TB')
ylabel('\chi/\mu')
set(gcf, 'PaperPosition', [0 0 1 2.5] * DrawParams.figW);
print(gcf, [Names.analysisDir, 'MFT_model.png'], '-dpng', ['-r',num2str(DrawParams.rez)], '-opengl')
close(gcf)

% time
dx        = SimParams.dx;
dt        = SimParams.dt;
nSteps    = round(SimParams.totTime/dt); % total number of simulation time steps
nOutSteps = round(SimParams.outDt/dt); % number of simulation time steps that produces an output
nOut      = floor(nSteps/nOutSteps); % number of outputs
tt        = 1; % output time index

% channel set up
L = SimParams.L; % number of spatial index
x = SimParams.x; % um, position along channel
if SimParams.dim == 2
    r = SimParams.r; % um, radial position in 2D
end

% cell
Drho = SimParams.Drho; % cell spatial diffusion coefficient
DL   = SimParams.DL;
% centrifuge parameters
if SimParams.centrifuge == 1
    r0 = SimParams.r0; % distance in um from the beginning of the channel to the center of the SimParams.centrifuge, 7 cm
    w0 = SimParams.w0; % spinning speed, 3000 rpm
    sedCoeff = SimParams.sedCoeff; % sedimentation coefficient of E. coli
    ux = ones(1, DL);
else
    ux = exp(-(1:2*DL).^2/(DL)^2);
end
ux = ux./sum(ux); % initial cell probability
initial_total_density = (CW.P*CW.dCWBias*ux)*SimParams.NCells/SimParams.OD2Cell; % OD
if SimParams.secondPop.phi > 0
    u = zeros(CW.ncwb, L, 2); % OD, cell density with dimension CW.CWBias * x * [pop1 pop2]
    u(:,2:(2*DL+1),2) = initial_total_density * SimParams.secondPop.phi; % update cell density of the second population
    aA = [aA aA*SimParams.secondPop.aARatio]';
    if SimParams.oxy > 0
        aO = [aO aO*SimParams.secondPop.aORatio]';
    end
    if SimParams.ser > 0
        aS = [aS aS*SimParams.secondPop.aSRatio]';
    end
    if SimParams.nut > 0
        aN = [aN aN*SimParams.secondPop.aNRatio]';
    end
    if SimParams.glu > 0
        sG = [sG sG*SimParams.secondPop.sGRatio]';
        cG = [cG cG*SimParams.secondPop.cGRatio]';
    end
else
    u = zeros(CW.ncwb, L); % OD, cell density with dimension CW.CWBias * x * [1]
end

% predict theory and update cell density of the first population
if SimParams.concType == 3
    WaveResults = [];
    WaveResults.z = dx + (0:dx:5000);
    WaveResults.c = SimParams.aA * sum(sum(initial_total_density)) * dx / SimParams.asp;
    WaveResults.tbDistAll = CW.P;
    WaveResults.dynDen = initial_total_density;
    WaveResults.tbMaxAtLastEps = 1;
    WavePrediction = shapeprediction(WaveResults, SimParams, CW, MFTmodel);
    predSize = numel(WaveResults.z);
    
    % use prediction to update cell density
    if SimParams.secondPop.phi > 0 % only use theory when there is a single population
        error('Cannot have SimParams.concType == 3 while SimParams.secondPop.phi > 0')
    end
    for i = 1:CW.ncwb
        [~, indZMax] = min(abs(WavePrediction.zMaxPred(i) - WaveResults.z));
        u(i,indZMax+1) = CW.P(i)*CW.dCWBias*SimParams.NCells/SimParams.OD2Cell; % all cells of phenotype i are concentrated at the z-block indZMax+1
    end
elseif SimParams.concType == 4
    if SimParams.secondPop.phi > 0 % only use scaling when there is a single population
        error('Cannot have SimParams.concType == 4 while SimParams.secondPop.phi > 0')
    end
    DL = round(SimParams.xScale/SimParams.dx/5);
    u(:,2:DL+1) = CW.P*CW.dCWBias*ones(1,DL)/DL*SimParams.NCells/SimParams.OD2Cell;
else
    u(:,2:(2*DL+1),1) = initial_total_density * (1-SimParams.secondPop.phi);
end

% finish cell density initialization with ghost cells
u(:,1,:) = u(:,2,:);

figure('visible','off')
plot(x, sum(u(:,:,1)), 'color',[0. 0.5 0.])
legend('cells')
if SimParams.secondPop.phi > 0
    hold on
    plot(x, sum(u(:,:,2)), 'r')
    legend('population 1', 'population 2')
end
xlabel('x (\mum)')
ylabel('cell density (OD)')
set(gcf, 'PaperPosition', [0 0 1 0.7] * DrawParams.figW);
print(gcf, [Names.analysisDir, 'initial_cells.png'], '-dpng', ['-r',num2str(DrawParams.rez)], '-opengl')
close(gcf)

% aspartate
DA   = SimParams.DA; % um^2/s, aspartate diffusion coefficient, A method for the determination of diffusion coefficients for small molecules in aqueous solution, Analytical Biochemistry, Volume 166, Issue 2, 1 November 1987, Pages 335-341
aA   = SimParams.aA; % uM/OD/s, asp consumption rate when cell density is low
KA   = SimParams.KA; % uM, Michaelis-Menten constant in asp consumption
lAO  = SimParams.lAO; % fraction of asp consumption without oxygen
KAO  = SimParams.KAO; % uM, scaling constant in oxygen-dependence of asp consumption
HAO  = SimParams.HAO; % steepness in oxygen-dependence of asp consumption
switch SimParams.concType
    case {0, 2}
        A = SimParams.asp*1./(1+(SimParams.V./x).^SimParams.H); % uM, initial sigmoidal aspartate concentration
    case {1, 4}
        A = ones(1,L)*SimParams.asp;
    case 3
        A = ones(1,L)*SimParams.asp;
        A(2:predSize+1) = WavePrediction.aspPred;
        A(1) = A(2);
end

figure('visible','off')
plot(x, A)
xlabel('x (\mum)')
ylabel('Asp (\muM)')
set(gcf, 'PaperPosition', [0 0 1 0.7] * DrawParams.figW);
print(gcf, [Names.analysisDir, 'initial_asp.png'], '-dpng', ['-r',num2str(DrawParams.rez)], '-opengl')
close(gcf)

% oxygen
if SimParams.oxy > 0
    DO = SimParams.DO; % um^2/s, oxygen diffusion coefficient, http://compost.css.cornell.edu/oxygen/oxygen.diff.water.html
    aO = SimParams.aO; % uM/OD/s, maximal oxygen consumption rate
    KO = SimParams.KO; % uM, Michaelis-Menten constant in oxygen consumption
    kO = SimParams.kO; % 1/s, oxygen transfer rate from PDMS to liquid interface
    switch SimParams.concType
        case {0, 2}
            O = SimParams.oxy*1./(1+(SimParams.V./x).^SimParams.H); % uM, initial sigmoidal oxygen concentration
        case {1, 4}
            O = ones(1,L)*SimParams.oxy;
        case 3
        	O = ones(1,L)*SimParams.oxy;
            O(2:predSize+1) = WavePrediction.oxyPredExact;
    end
    
    figure('visible','off')
    plot(x, O)
    xlabel('x (\mum)')
    ylabel('Oxy (\muM)')
    set(gcf, 'PaperPosition', [0 0 1 0.7] * DrawParams.figW);
    print(gcf, [Names.analysisDir, 'initial_oxy.png'], '-dpng', ['-r',num2str(DrawParams.rez)], '-opengl')
    close(gcf)
else
    O = zeros(1,L);
end

% serine
if SimParams.ser > 0
    DS = SimParams.DS; % um^2/s, serine diffusion coefficient
    aS = SimParams.aS; % uM/OD/s, maximal serine consumption rate
    KS = SimParams.KS; % uM, Michaelis-Menten constant in serine consumption
    switch SimParams.concType
        case {0, 2}
            S = SimParams.ser*1./(1+(SimParams.V./x).^SimParams.H); % uM, initial sigmoidal serine concentration
        case {1, 4}
            S = ones(1,L)*SimParams.ser;
        case 3
            error('Cannot have SimParams.concType == 3 while SimParams.ser > 0')
    end
    
    figure('visible','off')
    plot(x, S)
    xlabel('x (\mum)')
    ylabel('Ser (\muM)')
    set(gcf, 'PaperPosition', [0 0 1 0.7] * DrawParams.figW);
    print(gcf, [Names.analysisDir, 'initial_ser.png'], '-dpng', ['-r',num2str(DrawParams.rez)], '-opengl')
    close(gcf)
else
    S = zeros(1,L);
end

% nutrient
if SimParams.nut > 0
    DN = SimParams.DN; % um^2/s, nutrient diffusion coefficient
    aN = SimParams.aN; % uM/OD/s, nutrient consumption rate
    KN = SimParams.KN; % uM, Michaelis-Menten constant in nutrient consumption
    switch SimParams.concType
        case {0, 2}
            N = SimParams.nut*1./(1+(SimParams.V./x).^SimParams.H); % mg/mL, initial sigmoidal nutrient concentration
        case {1, 4}
            N = ones(1,L)*SimParams.nut;
        case 3
            error('Cannot have SimParams.concType == 3 while SimParams.nut > 0')
    end
    
    figure('visible','off')
    plot(x, N)
    xlabel('x (\mum)')
    ylabel('Nutrient (mg/mL)')
    set(gcf, 'PaperPosition', [0 0 1 0.7] * DrawParams.figW);
    print(gcf, [Names.analysisDir, 'initial_nutrient.png'], '-dpng', ['-r',num2str(DrawParams.rez)], '-opengl')
    close(gcf)
else
    N = zeros(1,L);
end

% glutamate
if SimParams.glu > 0
    DG = SimParams.DG; % um^2/s, glutamate diffusion coefficient
    sG = SimParams.sG; % uM/OD/s, glutamate secretion rate
    dG = SimParams.dG; % uM/OD/s, glutamate degradation rate
    cG = SimParams.cG; % uM/OD/s, glutamate consumption rate
end
G = zeros(1,L);

    function fS = receptorSensing(L, Ki, Ka) % convert to perceived signal
        fS2 = L/Ki; % log((1+S/Ki/100));
        fS1 = log((1+10000*L/Ki)); % log(fS2);
        fS0 = log((1+L/Ki)./(1+L/Ka));
        switch SimParams.recType
            case 0
                fS = fS0;
            case 1
                fS = fS1;
            case 2
                fS = fS2;
        end
        if tp
            figure(tp2)
            showInds = 1:50;
            subplot(3,1,1)
            plot(1000*x(showInds), fS0(:,showInds), '.')
            xlabel('x (mm)')
            ylabel('\phi_N_L')
            subplot(3,1,2)
            plot(1000*x(showInds), fS1(:,showInds), '.')
            xlabel('x (mm)')
            ylabel('\phi_L_G')
            subplot(3,1,3)
            plot(1000*x(showInds), fS2(:,showInds), '.')
            xlabel('x (mm)')
            ylabel('\phi_L_I')
        end
    end

    function [ui, Ai, Oi, Si, Ni, Gi] = integrateStep(u1, A1, O1, S1, N1, G1) % atomic integration step
        % initialize
        ui = zeros(size(u1));
        Ai = zeros(size(A1));
        Oi = zeros(size(O1));
        Si = zeros(size(S1));
        Ni = zeros(size(N1));
        Gi = zeros(size(G1));
        
        % cell density dynamics
        % calculate perceived signal
        f = receptorSensing(A1, SimParams.KiA, SimParams.KaA);
        if SimParams.oxy > 0
           f = f + SimParams.chiO * receptorSensing(O1, SimParams.KiO, SimParams.KaO);
        end
        if SimParams.ser > 0
           f = f + SimParams.chiS * receptorSensing(S1, SimParams.KiS, SimParams.KaS);
        end
        if SimParams.glu > 0
           f = f + SimParams.chiG * receptorSensing(G1, SimParams.KiG, SimParams.KaG);
        end
        % calculate chemotactic drift
        g = zeros(size(u1)); % initialize drift term as temporary variable
        g(:,2:end-1,1) = MFTmodel.chi * (f(3:end) - f(1:end-2))/(2*dx); % chemotactic drift
        if SimParams.secondPop.phi > 0 % modifies drift term of the second population
            g(:,2:end-1,2) = g(:,2:end-1,1) * SimParams.secondPop.speedRatio^2; % chi ~ v^2
        end
        % calculate centrifugal drift
        if SimParams.centrifuge
            g(:,2:end-1,:) = g(:,2:end-1,:) - sedCoeff * w0^2 * (r0 - x(2:end-1)); % will subtract the second dimension copied into other dimensions
        end
        % summarize total drift
        g = g .* u1;
        g(:,1,:) = -g(:,2,:); % set flux at ghost cell and the first cell opposite to net 0 to create channel beginning
        g(:,end,:) = -g(:,end-1,:); % set flux at ghost cell and the last cell opposite to net 0 to create channel end
        ui(:,2:end-1,:) = Drho(:,2:end-1) .* (u1(:,1:end-2,:) + u1(:,3:end,:) - 2*u1(:,2:end-1,:))/dx^2 ... % diffusion coefficient will copy into the third dimension
                          - (g(:,3:end,:) - g(:,1:end-2,:))/(2*dx); % diffusion term
        % modify divergence in 2D
        if SimParams.dim == 2 % radius will copy from the second dimension into other dimensions
            ui(:,2:end-1,:) = ui(:,2:end-1,:) ...
                              + 1./r(2:end-1) .* Drho(:,2:end-1) .* (u1(:,3:end,:) - u1(:,1:end-2,:))/(2*dx) ...
                              - 1./r(2:end-1) .* g(:,2:end-1,:);
        end
        % add cell growth
        if SimParams.growthRate > 0 && SimParams.nut > 0 % growth is proportional to nutrient
            if SimParams.secondPop.phi > 0
                mutation = cat(3, SimParams.PNew * u1(:,2:end-1,1), SimParams.PNew * u1(:,2:end-1,2));
            else
                mutation = SimParams.PNew * u1(:,2:end-1);
            end
            ui(:,2:end-1,:) = ui(:,2:end-1,:) ...
                              + SimParams.growthRate * N1(:,2:end-1) .* ((SimParams.inherit * u1(:,2:end-1,:)) ...
                              + ((1 - SimParams.inherit) * (2 * mutation - u1(:,2:end-1,:))));
        end
        ui(:,1,:) = ui(:,2,:); % set density equal at ghost cell & the first cell
        ui(:,end,:) = ui(:,end-1,:); % set density equal at ghost cell & the last cell

        % aspartate concentration dynamics
        maximalConsumption = squeeze(sum(u1)) * aA;
%         GammaA = maximalConsumption(:)' ./ (1+KA./A1) .* ((KAO./(1+((SimParams.oxy-O1)/ktAO).^5))+1-KAO); % uM/s, aspartate consumption term
%         GammaA = maximalConsumption(:)' ./ (1+KA./A1) .* (lAO + (1-lAO)*1./(1+(KAO./O1).^HAO)); % uM/s, aspartate consumption term
        GammaA = maximalConsumption(:)' ./ (1+KA./A1) .* (lAO + (1-lAO)*(O1/KAO).^HAO); % uM/s, aspartate consumption term
        Ai(2:end-1) = DA * (A1(1:end-2) + A1(3:end) - 2*A1(2:end-1))/dx^2 - GammaA(2:end-1);
        if SimParams.dim == 2
                Ai(2:end-1) = Ai(2:end-1) + 1./r(2:end-1) * DA .* (A1(3:end) - A1(1:end-2))/(2*dx);
        end
        Ai(1) = Ai(2);
        Ai(end) = Ai(end-1);

        % oxygen concentration dynamics
        if SimParams.oxy > 0
            maximalConsumption = squeeze(sum(u1)) * aO;
            GammaO = maximalConsumption(:)' ./ (1+KO./O1); % uM/s, oxygen consumption term
            Oi(2:end-1) = DO * (O1(1:end-2) + O1(3:end) - 2*O1(2:end-1))/dx^2 - GammaO(2:end-1) + kO * (SimParams.oxy - O1(2:end-1));
            if SimParams.dim == 2
                Oi(2:end-1) = Oi(2:end-1) + 1./r(2:end-1) * DO .* (O1(3:end) - O1(1:end-2))/(2*dx);
            end
            Oi(1) = Oi(2);
            Oi(end) = Oi(end-1);
        end
        
        % serine concentration dynamics
        if SimParams.ser > 0
            maximalConsumption = squeeze(sum(u1)) * aS;
            GammaS = maximalConsumption(:)' ./ (1+KS./S1); % uM/s, serine consumption term
            Si(2:end-1) = DS * (S1(1:end-2) + S1(3:end) - 2*S1(2:end-1))/dx^2 - GammaS(2:end-1);
            if SimParams.dim == 2
                Si(2:end-1) = Si(2:end-1) + 1./r(2:end-1) * DS .* (S1(3:end) - S1(1:end-2))/(2*dx);
            end
            Si(1) = Si(2);
            Si(end) = Si(end-1);
        end
        
        % nutrient concentration dynamics\
        if SimParams.nut > 0
            maximalConsumption = squeeze(sum(u1)) * aN;
            GammaN = maximalConsumption(:)' ./ (1+KN./N1); % mg/mL/s, constant nutrient consumption rate
            Ni(2:end-1) = DN * (N1(1:end-2) + N1(3:end) - 2*N1(2:end-1))/dx^2 - GammaN(2:end-1);
            if SimParams.dim == 2
                Ni(2:end-1) = Ni(2:end-1) + 1./r(2:end-1) * DN .* (N1(3:end) - N1(1:end-2))/(2*dx);
            end
            Ni(1) = Ni(2);
            Ni(end) = Ni(end-1);
        end
        
        % glutamate concentration dynamics
        if SimParams.glu > 0
            density = squeeze(sum(u1));
            secretionG = density * sG;
            secretionG = secretionG(:)'; % uM/s, glutamate secretion term
            degradationG = dG * G1; % uM/s, glutamate degradation term
            consumptionG = density * cG;
            consumptionG = consumptionG(:)' .* G1; % uM/s, glutamate consumption term
            Gi(2:end-1) = DG * (G1(1:end-2) + G1(3:end) - 2*G1(2:end-1))/dx^2 + secretionG(2:end-1) - consumptionG(2:end-1) - degradationG(2:end-1);
            if SimParams.dim == 2
                Gi(2:end-1) = Gi(2:end-1) + 1./r(2:end-1) * DG .* (G1(3:end) - G1(1:end-2))/(2*dx);
            end
            Gi(1) = Gi(2);
            Gi(end) = Gi(end-1);
        end
    end

%% simulation steps
tp = false; % test plot option for debug
if tp
    tp1 = figure;
    tp2 = figure;
end

% initialize temporal profiles
if SimParams.secondPop.phi > 0
    ut = zeros(CW.ncwb, L, nOut, 2);
else
    ut = zeros(CW.ncwb, L, nOut);
end
At = zeros(L, nOut);
Ot = zeros(L, nOut);
St = zeros(L, nOut);
Nt = zeros(L, nOut);
Gt = zeros(L, nOut);

figure%('visible','off')

for t = 1:nSteps % Runge-Kutta RK4, https://en.wikipedia.org/wiki/Runge%E2%80%93Kutta_methods
    % normalize cell counts in the wave
    if SimParams.cellNormalize && t == round(SimParams.normalizeTime / dt)
        xIndStart = round(SimParams.normalizePosition / SimParams.dx);
        currentNCells = sum(sum(u(:,xIndStart:end,:))) * SimParams.OD2Cell;
        scalingFactor = SimParams.normalizeNCells / currentNCells;
        disp(['Scaling up cell density at t = ', num2str(t*dt), ' s by a factor of ', num2str(scalingFactor), ' ...'])
        u = u * scalingFactor;
    end
    
    % test plot
    if tp
        figure(tp1)
        showInds = 1:50;
        subplot(2,1,1)
        plot(x(showInds), sum(sum(u(:,showInds,:))), '.')
        xlabel('x (mum)')
        ylabel('u (OD)')
        subplot(2,1,2)
        plot(x(showInds), A(showInds), '.')
        xlabel('x (\mum)')
        ylabel('A (uM)')
    end
    
    % Intermediate variable 1
    u1 = u;
    u1(u1<0) = 0;
    A1 = A;
    A1(A1<0) = 0;
    O1 = O;
    O1(O1<0) = 0;
    S1 = S;
    S1(S1<0) = 0;
    N1 = N;
    N1(N1<0) = 0;
    G1 = G;
    G1(G1<0) = 0;
    [ua, Aa, Oa, Sa, Na, Ga] = integrateStep(u1, A1, O1, S1, N1, G1);
    
    % Intermediate variable 2
    u1 = u + ua*dt/2;
    u1(u1<0) = 0;
    A1 = A + Aa*dt/2;
    A1(A1<0) = 0;
    O1 = O + Oa*dt/2;
    O1(O1<0) = 0;
    S1 = S + Sa*dt/2;
    S1(S1<0) = 0;
    N1 = N + Na*dt/2;
    N1(N1<0) = 0;
    G1 = G + Ga*dt/2;
    G1(G1<0) = 0;
    [ub, Ab, Ob, Sb, Nb, Gb] = integrateStep(u1, A1, O1, S1, N1, G1);
    
    % Intermediate variable 3
    u1 = u + ub*dt/2;
    u1(u1<0) = 0;
    A1 = A + Ab*dt/2;
    A1(A1<0) = 0;
    O1 = O + Ob*dt/2;
    O1(O1<0) = 0;
    S1 = S + Sb*dt/2;
    S1(S1<0) = 0;
    N1 = N + Nb*dt/2;
    N1(N1<0) = 0;
    G1 = G + Gb*dt/2;
    G1(G1<0) = 0;
    [uc, Ac, Oc, Sc, Nc, Gc] = integrateStep(u1, A1, O1, S1, N1, G1);
    
    % Intermediate variable 4
    u1 = u + uc*dt;
    u1(u1<0) = 0;
    A1 = A + Ac*dt;
    A1(A1<0) = 0;
    O1 = O + Oc*dt;
    O1(O1<0) = 0;
    S1 = S + Sc*dt;
    S1(S1<0) = 0;
    N1 = N + Nc*dt;
    N1(N1<0) = 0;
    G1 = G + Gc*dt;
    G1(G1<0) = 0;
    [ud, Ad, Od, Sd, Nd, Gd] = integrateStep(u1, A1, O1, S1, N1, G1);
    
    % Combine
    u = u + (ua + 2*ub + 2*uc + ud)*dt/6;
    A = A + (Aa + 2*Ab + 2*Ac + Ad)*dt/6;
    O = O + (Oa + 2*Ob + 2*Oc + Od)*dt/6;
    S = S + (Sa + 2*Sb + 2*Sc + Sd)*dt/6;
    N = N + (Na + 2*Nb + 2*Nc + Nd)*dt/6;
    G = G + (Ga + 2*Gb + 2*Gc + Gd)*dt/6;
    u(:,1,:) = u(:,2,:);
    u(:,end,:) = u(:,end-1,:);
    u(u<0) = 0;
    A(1) = A(2);
    A(end) = A(end-1);
    A(A<0) = 0;
    O(1) = O(2);
    O(end) = O(end-1);
    O(O<0) = 0;
    S(1) = S(2);
    S(end) = S(end-1);
    S(S<0) = 0;
    N(1) = N(2);
    N(end) = N(end-1);
    N(N<0) = 0;
    G(1) = G(2);
    G(end) = G(end-1);
    G(G<0) = 0;

    % find peak
    indx = find(diff(sign( diff(sum(sum(u,3))) )) == -2)+1; % find where first derivative = 0 and second derivative < 0
    if isempty(indx)
        indx = 1; 
    end
    
    % output profiles at regular time intervals
    if (mod(t,nOutSteps)==0)
        disp([num2str(t),' out of ', num2str(nSteps)])

        % save
        ut(:,:,tt,:) = u;
        uSum = squeeze(sum(u));
        At(:,tt)   = A;
        if SimParams.oxy > 0
            Ot(:,tt) = O;
        end
        if SimParams.ser > 0
            St(:,tt) = S;
        end
        if SimParams.dim == 2
            Nt(:,tt) = N;
        end
        if SimParams.glu > 0
            Gt(:,tt) = G;
            subplotTotal = 3;
        else
            subplotTotal = 2;
        end
        tt = tt+1;
        
        subplot(subplotTotal,1,1);
        hold on
        if SimParams.secondPop.phi > 0
            plot(x, uSum(:,1), 'color', [0. 0.5 0.])
            plot(x, uSum(:,2), 'r')
            legend('population1', 'population2')
        else
            plot(x, uSum, 'color', [0. 0.5 0.])
            legend('cells')
        end
        axis tight
        xlabel('x (\mum)')
        ylabel('\rho (OD)')
        
        subplot(subplotTotal,1,2)
        hold on
        plot(x,A/SimParams.asp,'b')
        lgs = {'asp'};
        if SimParams.oxy > 0
            plot(x,O/SimParams.oxy,'k')
            lgs = [lgs 'oxy'];
        end
        if SimParams.ser > 0
            plot(x,S/SimParams.ser,'m')
            lgs = [lgs 'ser'];
        end
        if SimParams.dim == 2
            plot(x,N/SimParams.nutrient,'c')
            lgs = [lgs 'nutrient'];
        end
        legend(lgs)
        axis tight
        xlabel('x (\mum)')
        ylabel('normalized concentrations')
        
        if SimParams.glu > 0
            subplot(subplotTotal,1,3)
            hold on
            plot(x, G, 'Color', [0.8 0.3 0])
            legend('glu')
            ylim([0 2000])
            xlabel('x (\mum)')
            ylabel('glu (\muM)')
        end
    end
end

set(gcf, 'PaperPosition', [0 0 1 subplotTotal*0.3] * DrawParams.figW);
saveas(gcf, [Names.analysisDir, 'wave_profiles.fig'])
print(gcf, [Names.analysisDir, 'wave_profiles.png'], '-dpng', ['-r',num2str(DrawParams.rez)], '-opengl')
close(gcf)

%% save output
SimResults.finalProfile.cellDensity  = squeeze(u); 
SimResults.dynProfile  .cellDensity  = squeeze(ut); % CW.ncwb * x * SimParams.time (* 2) array of cell density in OD
SimResults.finalProfile.asp          = A;
SimResults.dynProfile  .asp          = At; % x * SimParams.time array of asp concentration in uM
SimResults.finalProfile.oxy          = O;
SimResults.dynProfile  .oxy          = Ot; % x * SimParams.time array of oxy concentration in uM
SimResults.finalProfile.ser          = S;
SimResults.dynProfile  .ser          = St; % x * SimParams.time array of ser concentration in uM
SimResults.finalProfile.nut          = N;
SimResults.dynProfile  .nut          = Nt; % x * SimParams.time array of oxy concentration in mg/mL
SimResults.finalProfile.glu          = G;
SimResults.dynProfile  .glu          = Gt; % x * SimParams.time array of glu concentration in uM
end
%% Main function. Wrapper to run multiple simulations in parallel
%This script has no imputs. It just runs simulations by calling simulate
%wave, a function which takes in the struct SimParams. Various parameters
%can be edited after calling parameters(). 

clear 
close all

%% I vs. Chi
%List of parameter values to try
%Values to scale cell2 I and Chi by.
cell2I = [1];
cell2Chi = [0.18, 0.42, 1, 2.37, 5.62];

tic()
parfor ii = 1:length(cell2I)
    for jj = 1:length(cell2Chi)    
        %Get parameters from file
        SimParams = parameters();
        SimParams.irate2 = SimParams.irate2 .* cell2I(ii);
        SimParams.Chi2 = SimParams.Chi2 .* cell2Chi(jj);
        SimParams.OutFolderName = './Outputs/OutputsIChi/';
	SimParams.PFU = 0;	
        
        %Run simulation
        simulateWave(SimParams)
    end
end
toc()


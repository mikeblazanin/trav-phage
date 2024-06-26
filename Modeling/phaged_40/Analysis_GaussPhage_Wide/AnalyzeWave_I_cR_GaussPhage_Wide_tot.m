%% Load Wave Data
clear all
close all

FileList = dir('../Outputs_WideGaussPhage/OutputsIcR_GaussPhage/*.mat');

%Store measurements (PoplationCell 1, Pop. Cell 2, relative Chi,
%relative i)
Measurements = zeros(length(FileList), 6); 

for ii = 1:length(FileList)
    File = ['../Outputs_WideGaussPhage/OutputsIcR_GaussPhage/', FileList(ii).name];
    load(File)

    %% Get wave at last timepoint
    rho_cell_end = rho_cell_store(end, :);
    rho_cell2_end = rho_cell2_store(end, :);

    % set(gca', 'yscale', 'log')

%     %% Get indices of the wave
%     %Find the peaks and start of wave 1
%     peaks = findpeaks(rho_cell_end);
%     peakInd = find(rho_cell_end == peaks(end));
%     WaveStartInd = max(find(rho_cell_end(1:peakInd) <= min(rho_cell_end(1:peakInd))));
% 
%     %Find peaks and start of wave 2
%     peaks2 = findpeaks(rho_cell2_end);
%     peakInd2 = find(rho_cell2_end == peaks2(end));
%     WaveStartInd2 = max(find(rho_cell2_end(1:peakInd2) <= min(rho_cell2_end(1:peakInd))));

%     figure();
%     plot(rho_cell_end)
%     hold on
%     xline(peakInd)
%     xline(WaveStartInd)
% 
%     plot(rho_cell2_end)
%     xline(peakInd2)
%     xline(WaveStartInd2)

    %% Sum cell density in the wave
    cell_population = sum(rho_cell_end(1:end));
    cell2_population = sum(rho_cell2_end(1:end));
    
    %% Get Parameters
    relativecR = SimParams.cR2 ./ SimParams.cR;
    RelativeI = SimParams.irate2 ./ SimParams.irate;
    
    %% Store measurements
    Measurements(ii, :) = [cell_population, cell2_population, relativecR, RelativeI, SimParams.cR, SimParams.irate];
end

%% Plot
relcR = Measurements(:, 3);
relI = Measurements(:, 4);
relPopulation = Measurements(:, 2) ./ Measurements(:, 1);


figure()
scatter(relcR, relI, 60, log(relPopulation), 'filled')
a = colorbar;
colormap jet
a.Label.String = 'log N1/N2';
a.Label.Rotation = -90;
xlabel("Relative cR")
ylabel("Relative I")
set(gca,'xscale', 'log', 'yscale', 'log')
saveas(gcf, './I_vs_cR_gPhage_Wide_tot.png')

%% Export Data
headers = ["Cell_population", "Cell2_population", "relativecR", "relativeI", "cR", "irate"];
textHeader = strjoin(headers, ',');
fid = fopen('./I_vs_cR_gPhage_Wide_tot.csv', 'w');
fprintf(fid,'%s\n',textHeader)
fclose(fid)

dlmwrite('./I_vs_cR_gPhage_Wide_tot.csv',Measurements,'-append');

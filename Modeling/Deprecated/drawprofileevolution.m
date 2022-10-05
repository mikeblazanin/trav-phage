function drawprofileevolution(AnalysisParams, SimResults, SimParams, DrawParams, Names)
% This function draws temporal evolution of wave profiles.

% OUTPUT
% profile_evolution.fig
% profile_evolution.png

x = SimParams.x;
outDt = SimParams.outDt;
L = SimParams.L;

disp('Drawing profile evolution ...')

% plot normalized profile evolutions of all variables
close all
figure('visible', 'off')
subplot(2,2,1)
hold on
for t = 1:SimParams.nT
    cellDensity = squeeze(sum(SimResults.dynProfile.cellDensity(:,:,t,:)));
	maxCellDensity = max(cellDensity(:));
    maxGlutamate = max(SimResults.dynProfile.glu(:));
    
    plot3(x, t*outDt*ones(1,L), SimResults.dynProfile.asp(:,t)/SimParams.asp, 'b')
    lgs = {'asp'};
    
    if SimParams.oxy > 0
        plot3(x, t*outDt*ones(1,L), SimResults.dynProfile.oxy(:,t)/SimParams.oxy, 'k')
        lgs = [lgs 'oxy'];
    end
    
    if SimParams.glu > 0
        plot3(x, t*outDt*ones(1,L), SimResults.dynProfile.glu(:,t)/maxGlutamate, 'Color', [0.8 0.3 0])
        lgs = [lgs 'glu'];
    end
    
    if SimParams.ser > 0
        plot3(x, t*outDt*ones(1,L), SimResults.dynProfile.ser(:,t)/SimParams.ser, 'm')
        lgs = [lgs 'ser'];
    end
    
    if SimParams.dim == 2
        plot3(x, t*outDt*ones(1,L), SimResults.dynProfile.nut(:,t)/SimParams.nut, 'c')
        lgs = [lgs 'nut'];
    end
    
    if SimParams.secondPop.phi > 0
        plot3(x, t*outDt*ones(1,L), cellDensity(:,1)/maxCellDensity, 'color',[0. 0.5 0.])
    	plot3(x, t*outDt*ones(1,L), cellDensity(:,2)/maxCellDensity, 'r')
        lgs = [lgs 'population 1' 'population 2'];
    else
        plot3(x, t*outDt*ones(1,L), cellDensity/maxCellDensity, 'color',[0. 0.5 0.])
        lgs = [lgs 'cells'];
    end
    
    legend(lgs)
end

xlabel('x (\mum)')
ylabel('t (s)')
zlabel('normalized quantities')
view(-20,30)

% plot colormap of cell density
twoCellDensity = squeeze(sum(SimResults.dynProfile.cellDensity, 1));
densityColor = zeros(size(twoCellDensity,1), size(twoCellDensity,2), 3);
if ndims(twoCellDensity) == 2
    densityColor(:,:,2) = twoCellDensity/max(twoCellDensity(:));
else
    densityColor(:,:,2) = twoCellDensity(:,:,1)/max(max(twoCellDensity(:,:,1))) .* (twoCellDensity(:,:,1) > twoCellDensity(:,:,2));
    densityColor(:,:,1) = twoCellDensity(:,:,2)/max(max(twoCellDensity(:,:,2))) .* (twoCellDensity(:,:,1) <= twoCellDensity(:,:,2));
end
densityColor = flip(densityColor);    

subplot(2,2,3)
image([SimParams.time(1), SimParams.time(end)], [SimParams.x(end), SimParams.x(1)], densityColor)
set(gca,'YDir','normal')
hold on
plot(SimParams.outDt * size(densityColor,2)/2 * [1 1], SimParams.dx * [0 size(densityColor,1)], 'w--', 'LineWidth', DrawParams.lw)
xlabel('t (s)')
ylabel('x (\mum)')

axes('Position',[0.13 0.32 .15 .13])
hold on
if ndims(twoCellDensity) == 2
    plot(twoCellDensity(:, ceil(size(densityColor,2)/2)), 'color',[0. 0.5 0.], 'LineWidth', DrawParams.lw)
else
    plot(twoCellDensity(:, ceil(size(densityColor,2)/2), 1), 'color',[0. 0.5 0.], 'LineWidth', DrawParams.lw)
    plot(twoCellDensity(:, ceil(size(densityColor,2)/2), 2), 'r', 'LineWidth', DrawParams.lw)
end
set(gca, 'XTickLabel', {})
set(gca, 'YTickLabel', {})

% plot wave speed evolution
totalCellDensity = squeeze(sum(sum(SimResults.dynProfile.cellDensity, 4), 1));
c = nan(2, SimParams.nT); % array of wave speeds
for t = AnalysisParams.calculatePeakSpacing+1:SimParams.nT
    dsdI = diff(sign( diff(totalCellDensity(:,t-AnalysisParams.calculatePeakSpacing)) ));
    dsdF = diff(sign( diff(totalCellDensity(:,t)) ));
    denMaxIndI = find(dsdI == -2) + 1; % find where first derivative = 0 and second derivative < 0
    denMaxIndF = find(dsdF == -2) + 1; % find where first derivative = 0 and second derivative < 0
    if isempty(denMaxIndI)
        denMaxIndI1 = 1;
        denMaxIndI2 = 1;
    elseif numel(denMaxIndI) == 1
        denMaxIndI1 = denMaxIndI;
        denMaxIndI2 = denMaxIndI;
    else
        denMaxIndI1 = denMaxIndI(end);
        denMaxIndI2 = denMaxIndI(end-1);
    end
    if isempty(denMaxIndF)
        denMaxIndF1 = 1;
        denMaxIndF2 = 1;
    elseif numel(denMaxIndF) == 1
        denMaxIndF1 = denMaxIndF;
        denMaxIndF2 = denMaxIndF;
    else
        denMaxIndF1 = denMaxIndF(end);
        denMaxIndF2 = denMaxIndF(end-1);
    end
    
    % find wave speeds
    c(1,t) = (denMaxIndF1 - denMaxIndI1)/AnalysisParams.calculatePeakSpacing * SimParams.dx / SimParams.outDt;
    c(2,t) = (denMaxIndF2 - denMaxIndI2)/AnalysisParams.calculatePeakSpacing * SimParams.dx / SimParams.outDt;
end

subplot(2,2,2)
hold on
if ndims(twoCellDensity) ~= 2
    plot(SimParams.time, c(1,:), 'r', 'LineWidth', DrawParams.lw)
    plot(SimParams.time, c(2,:), 'color',[0. 0.5 0.], 'LineWidth', DrawParams.lw)
else
    plot(SimParams.time, c(1,:), 'color',[0. 0.5 0.], 'LineWidth', DrawParams.lw)
end
xlim([0 SimParams.time(end)])
ylim([0 DrawParams.cMax])
xlabel('t (s)')
ylabel('c (\mum/s)')
subplot(2,2,4)
dlogc1dt5 = -[0 diff(smooth(c(1,:),5)')]./c(1,:)/SimParams.outDt;
dlogc1dt10 = -[0 diff(smooth(c(1,:),10)')]./c(1,:)/SimParams.outDt;
dlogc1dt20 = -[0 diff(smooth(c(1,:),20)')]./c(1,:)/SimParams.outDt;
hold on
plot(SimParams.time, dlogc1dt5, ':', 'color',[0. 0.5 0.], 'LineWidth', DrawParams.lw)
plot(SimParams.time, dlogc1dt10, '--', 'color',[0. 0.5 0.], 'LineWidth', DrawParams.lw)
plot(SimParams.time, dlogc1dt20, '-', 'color',[0. 0.5 0.], 'LineWidth', DrawParams.lw)
if ndims(twoCellDensity) ~= 2
    dlogc2dt5 = -[0 diff(smooth(c(2,:),5)')]./c(2,:)/SimParams.outDt;
    dlogc2dt10 = -[0 diff(smooth(c(2,:),10)')]./c(2,:)/SimParams.outDt;
    dlogc2dt20 = -[0 diff(smooth(c(2,:),20)')]./c(2,:)/SimParams.outDt;
    plot(SimParams.time, dlogc2dt5, ':r', 'LineWidth', DrawParams.lw)
    plot(SimParams.time, dlogc2dt10, '--r', 'LineWidth', DrawParams.lw)
    plot(SimParams.time, dlogc2dt20, '-r', 'LineWidth', DrawParams.lw)
end
xlim([0 SimParams.time(end)])
ylim([0 1e-3])
legend('smooth 5', 'smooth 10', 'smooth 20')
xlabel('t (s)')
ylabel('dlog(c)/dt (1/s)')

set(gcf, 'PaperPosition', [0 0 1 0.6] * 2 * DrawParams.figW);
figName = [Names.analysisDir, 'profile_evolution'];
disp(['Saving figure ', figName, '...'])
saveas(gcf, [figName, '.fig'])
print(gcf, [figName, '.png'], '-dpng', ['-r',num2str(DrawParams.rez)], '-opengl')
close(gcf)
end
%% scan all external Asp concentrations
Asps = 10:10:200; % uM , total Asp concentration

% initialization
channel_length = 30000; % um
Tot_Time = 30*60; % s, total simulation time
out_dt = 30; % s, the output time interval
dx = 50; % um
dt = dx^2/5000; % s, relation with dx to ensure stability
OD1 = 8.4*10^8; % cells/ml in OD1
OD2Cell = dx*600*14*OD1/10000; % each block has volume dx*600um*14um
conc_type = 0; % initial concentration of Asp and Oxy, 0 = Hill-shaped, 1 = flat
ptb_type = 3; % initial tumble bias distribution, 0 = Bell-shaped, 1 = uniform, 2 = custom distribution, 3 = Adam's WT distribution
rec_type = 0; % receptor nonlinearity, 0 = nonlinear, 1 = log, 2 = linear
V = 800; % initial Asp profile length scale
H = 2.5; % initial Asp profile steepness
centrifuge = 0; % whether to include centrifugal force

% initial concentration of Asp and Oxy
switch conc_type
    case 0
        c_str = 'Hill';
    case 1
        c_str = 'Flat';
end

% initial tumble bias distribution
switch ptb_type
    case 0
        p_str = 'Bell';
    case 1
        p_str = 'Unif';
    case 2
        p_str = 'Cstm';
    case 3
        p_str = 'Adam';
end

% receptor nonlinearity
switch rec_type
    case 0
        r_str = 'NL';
    case 1
        r_str = 'LG';
    case 2
        r_str = 'LI';
end

% initialize cell model
[dCWBias, CWBias, P] = CWInitialize(ptb_type);
MFTmodel = MFTmodelInitialize(CWBias, ptb_type);

tb_dist = nan(numel(CWBias), numel(Asps));

for a = 20
    % names
    Asp = Asps(a);
    init_cond_dir = ['Asp_Oxy_', c_str, 'Conc_', p_str, 'TB_', r_str, 'rec\'];
    case_pref = ['wave_simulation_',num2str(Asp),'_',num2str(Tot_Time),...
        '_',num2str(out_dt),'_',num2str(channel_length),'_',num2str(dt),'_',num2str(dx)];
    sub_case_dir = '';
    case_dir = [case_pref, sub_case_dir, '\'];
    dropbox_dir = 'C:\Users\jl2345\Dropbox (emonetlab)\';
    data_dir = [dropbox_dir, 'users\junjiajia_long\data\SCRATCH\wave\', init_cond_dir, case_dir];
    analysis_dir = [dropbox_dir, 'users\junjiajia_long\analysis\wave\', init_cond_dir, case_dir];
    movie_dir = [analysis_dir, 'pheno_profile\'];
    sim_name = [data_dir, case_pref, '.mat'];
    analysis_name = [analysis_dir, case_pref, '_analysis.mat'];
    
    if exist(sim_name, 'file')
        % load
        disp(['Loading file ', sim_name, '...'])
        load(sim_name)
    else
        % run simulation
        disp(['Running simulation ', sim_name, '...'])
        [wavespeed_late, final_profile, dyn_profile, simParams]...
            = wave_simulation(Asp, Tot_Time, out_dt,...
            channel_length, dt, dx, conc_type, dCWBias, CWBias, P, MFTmodel,...
            V, H, centrifuge, rec_type);
        % save
        if ~exist(data_dir, 'dir')
            mkdir(data_dir)
        end
        disp(['Saving to file ', sim_name, '...'])
        save(sim_name, 'Asp', 'Tot_Time', 'out_dt', 'channel_length', 'dt', 'dx',...
            'MFTmodel', 'P', 'CWBias', 'dCWBias', 'simParams', 'conc_type', 'ptb_type', 'rec_type',...
            'wavespeed_late', 'final_profile', 'dyn_profile', '-v7.3')
    end
    
    % determine wave packet range
    cell_density = squeeze(sum(dyn_profile.cell_density(:, :, end))); % cell density in OD totaling all phenotypes
    dsd = diff(sign( diff(cell_density) ));
    den_min_ind = find(dsd == 2) + 1; % find where first derivative = 0 and second derivative > 0
    if isempty(den_min_ind)
        den_min_ind = 1;
    end
    den_min_ind = den_min_ind(1); % the first local minimum is real
    
    % find tb distribution within wave
    tb_dist(:,a) = sum(dyn_profile.cell_density(:, den_min_ind:end, end),2);
    tb_dist(:,a) = tb_dist(:,a)/sum(tb_dist(:,a))/dCWBias;
    
    % test figure
    figure
    subplot(2,1,1)
    plot(cell_density)
    subplot(2,1,2)
    hold on
    plot(CWBias, tb_dist(:,a))
    plot(CWBias, P)
    pause
end
%%
tb_dist_file = [dropbox_dir, 'users\junjiajia_long\analysis\wave\', init_cond_dir,'tb_dists.mat'];
save(tb_dist_file, 'tb_dist')
load(tb_dist_file)
%% plot tb
figure
hold on
for a = 1:numel(Asps)
    plot3(CWBias, Asps(a)*ones(size(CWBias)), tb_dist(:,a))
end
xlabel('TB')
ylabel('Asp (\muM)')
zlabel('PDF')
view(5,20)
c_tb = cumsum(tb_dist); % cumulative tb distribution
[xx, yy] = meshgrid(CWBias, Asps);
figure
hold on
contour(xx, yy, c_tb'*dCWBias, [0 0.7], 'k')
contour(xx, yy, c_tb'*dCWBias, [0 0.9], 'm')
contour(xx, yy, c_tb'*dCWBias, [0 0.95], 'c')
heatmap_plot = surf(xx, yy, c_tb'*dCWBias-1.5);
legend('70%', '90%', '95%', 'CDF')
set(heatmap_plot,'edgecolor','none')
view(2)
cb = colorbar;
caxis([-1.5 -0.5])
cb.TickLabels = {'0';'0.1';'0.2';'0.3';'0.4';'0.5';'0.6';'0.7';'0.8';'0.9';'1'};
title(cb,'CDF')
xlabel('TB')
ylabel('Asp (\muM)')
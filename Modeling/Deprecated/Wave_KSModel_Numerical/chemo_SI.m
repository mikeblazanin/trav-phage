% %%
% clear all;
% Aspind=log(30):(log(220)-log(30))/250:log(220);
% % Timeind=log(20):(log(55)-log(20))/100:log(55);

% tic;
% parfor i=1:length(Asp)
% clear all; close all;
% uiopen('D:\Dropbox\Yale\paper_wave_diversity\manuscript\Figures\final figure May2016\fig2\simvsexp\FigS7C_celldensity_v3.fig',1)

  Asp=[50 100 200];
   Time=[15 25 30];tic;

% clear all;
% Aspind=log(30):(log(220)-log(30))/50:log(220);
%  Asp=30:5:400;Time=12:43/(length(Asp)-1):55;
%  figure
  co=jet(length(Asp));
tic;

 dx=50/1000;
 Time0=0;
for j=3:3%length(Asp)
       [final_profile(j),wavespeed(j) dyn_profile(j)]=chemo_oxygen_fun_v4(Asp(j),Time(j)*60,Time0*60);
       cell_density=sum(final_profile(1).cell_density);
        x=0:dx:dx*(length(cell_density)-1);
        cell_density=sum(final_profile(j).cell_density);
         [maxden maxind]=max(cell_density(48:end));
%          subplot(2,1,1);
          
        plot(x-x(maxind+47),cell_density(1:end)*8.4*10^8,'Color',co(j,:));hold on;xlim([-1.5 1.5]);
%         subplot(2,1,2);
%         plot(peakpos,'Color',co(i,:));hold on;
end


%  
%   save D:\Asp_run.mat -v7.3
% % % % toc    
% % figure;
%  %% cell_density_profile_shifted averaged
% % figure;
% % dx=50;
% %   for j=1:length(Asp)
% %       s=size(dyn_profile(j).cell_density);
% %       cell_density=dyn_profile(j).cell_density;
% %       for i=5:s(3)
% %         cell_den=sum(cell_density(:,:,i)); 
% % %         x=0:4*dx:dx*(length(cell_den)-1);
% % %         cell_den=reshape(cell_den,4,length(cell_den)/4);cell_den=mean(cell_den);
% %         x=0:dx/1000:dx*(length(cell_den)-1)/1000;
% %         [maxden maxind]=max(cell_den(48:end));
% %         cell_den_shifted=cell_den(x-x(maxind+47)>-2&x-x(maxind+47)<4);
% %         cell_number(i-4)=sum(cell_den_shifted)*8.4*10^8*10^-12*600*14*dx;
% % %          plot(x-x(maxind+11),cell_den(1:end)*8.4*10^8);hold on;
% %       end
% %         cell_number_mean(j)=mean(cell_number);
% %         cell_number_std(j)=std(cell_number);
% %         clear cell_number;
% %   end
% % x=(1:length([41+maxind-40:41+maxind+40]))*dx-40*dx;
% % plot(x,cell_den_mean);
% %   xlim([-1.5 1.5]);
% 
%    %%
% % %     for j=1:3
% % %         cell_density=sum(final_profile(j).cell_density); x=0:4*dx:dx*(length(cell_density)-1);
% % %         cell_density=reshape(cell_density,4,length(cell_density)/4);cell_density=mean(cell_density);
% % %         [maxden maxind]=max(cell_density(12:end));
% % %         plot(x-x(maxind+11),cell_density(1:end)*8.4*10^8);hold on;
% % %   end
% % % 
% % %   xlim([-1.5 1.5]);
% %   
% %   %% total cell number by shifted average cell density profile;
% % dx=50;
% % for i=1:length(Asp)
% %     cell_number(i)=sum(cell_den_mean(:,i))*10^-12*600*14*dx;
% % end
% % loglog(Asp,cell_number(1:length(Asp))); hold on;
% % xlim([40 210]);
% % xlabel('Asp (\muM)');ylabel('cell number');
% % %% effective consumption rate
% %  for j=1:length(Asp)
% %       s=size(dyn_profile(j).cell_density);
% %       cell_density=dyn_profile(j).cell_density;
% %       for i=5:s(3)
% %         cell_den=sum(cell_density(:,:,i)); 
% % %         x=0:4*dx:dx*(length(cell_den)-1);
% % %         cell_den=reshape(cell_den,4,length(cell_den)/4);cell_den=mean(cell_den);
% %         [maxden maxind]=max(cell_den(42:end));
% %         cell_den_temp=cell_den(41+maxind-40:41+maxind+40);
% %         cell_number(i)=sum(cell_den_temp)*10^-12*600*14*dx;
% % %          plot(x-x(maxind+11),cell_den(1:end)*8.4*10^8);hold on;
% %       end
% % %         cell_den_mean(:,j)=mean(cell_den_temp')*8.4*10^8;
% % %         cell_den_std(:,j)=std(cell_den_temp')*8.4*10^8;
% % %         clear cell_den_temp;
% %         total_cell_num(j)=mean(cell_number);
% %         clear cell_number;
% %   end
% % consum_rate=Asp.*wavespeed./total_cell_num*10^-12*0.6*0.014;
% % plot(total_cell_num,consum_rate);
% %% tumble bias distribution
% % CWBias_exp=[0,0.0100000000000000,0.0200000000000000,0.0300000000000000,0.0400000000000000,0.0500000000000000,0.0600000000000000,0.0700000000000000,0.0800000000000000,0.0900000000000000,0.100000000000000,0.110000000000000,0.120000000000000,0.130000000000000,0.140000000000000,0.150000000000000,0.160000000000000,0.170000000000000,0.180000000000000,0.190000000000000,0.200000000000000,0.210000000000000,0.220000000000000,0.230000000000000,0.240000000000000,0.250000000000000,0.260000000000000,0.270000000000000,0.280000000000000,0.290000000000000,0.300000000000000,0.310000000000000,0.320000000000000,0.330000000000000,0.340000000000000,0.350000000000000,0.360000000000000,0.370000000000000,0.380000000000000,0.390000000000000,0.400000000000000,0.410000000000000,0.420000000000000,0.430000000000000,0.440000000000000,0.450000000000000,0.460000000000000,0.470000000000000,0.480000000000000,0.490000000000000,0.500000000000000,0.510000000000000,0.520000000000000,0.530000000000000,0.540000000000000,0.550000000000000,0.560000000000000,0.570000000000000,0.580000000000000,0.590000000000000,0.600000000000000,0.610000000000000,0.620000000000000,0.630000000000000,0.640000000000000,0.650000000000000,0.660000000000000,0.670000000000000,0.680000000000000,0.690000000000000,0.700000000000000];
% % P_exp=[0,0.00302823490236943,0.00498640060510877,0.00740654246062159,0.0129576724735464,0.0282364445536801,0.0342746634745129,0.0605832532805051,0.0797351468975508,0.105863629784042,0.144957679208212,0.222281447976504,0.281608928805501,0.373987712933409,0.512065751136908,0.691152325282671,0.919556834551468,1.19485305314468,1.51155357679430,1.88448868628358,2.26163802431774,2.67050119556754,3.04502589236555,3.42156794182481,3.78827243915744,4.11858437055020,4.26770099008759,4.46704480752076,4.50323811436766,4.50544966346908,4.42742856648661,4.31612752213579,4.11316467718262,3.94529362322380,3.65109481682462,3.23848249620281,2.98678093597204,2.69627429599337,2.36626462681502,2.10807472740901,1.95992530735338,1.71658120173055,1.53144380527692,1.43826939988074,1.27455605283792,1.14836389445003,1.01135546697596,0.907225653941910,0.780385038989034,0.685196631521439,0.645362573778435,0.620271119730271,0.579787157033761,0.525204939995039,0.523233085565830,0.466362640857505,0.445131198472857,0.399401650650870,0.399447904195710,0.371516048640824,0.380986255044554,0.393076006538824,0.375160618566594,0.357271202286408,0.341376560968248,0.322917046544370,0.291493170711647,0.311004430859554,0.290374229335908,0.272975724047314,0.162176873861474];
% % dCWBias=0.001;CWBias=0.001:dCWBias:CWBias_exp(end);
% % P=interp1(CWBias_exp,P_exp,CWBias,'spline');P=P./sum(P)./dCWBias;
% % % figure;
% % % plot(CWBias,P); hold on;
% % for i=1:length(Asp)
% %     cell_density=sum(final_profile(i).cell_density);   
% %     [maxden maxind]=max(cell_density(48:end));
% %     dyn_den_end=final_profile(i).cell_density(:,47+maxind-40:end);
% %     P_wave(:,i)=sum(dyn_den_end')./sum(sum(dyn_den_end'))./dCWBias;
% %     rel_wave(:,i)= P_wave(:,i)'./P;
% %     plot(CWBias,P_wave(:,i));
% %     hold on;
% % end
% %%  
% % figure;
% %   for j=1:3
% %         cell_density=sum(final_profile(j).cell_density); x=0:dx:dx*(length(cell_density)-1);
% % %         cell_density=reshape(cell_density,4,length(cell_density)/4);cell_density=mean(cell_density);
% %         [maxden maxind]=max(cell_density(48:end));
% %         plot(x-x(maxind+47),cell_density(1:end)./maxden);hold on;
% %   end
% % 
% %   xlim([-1.5 1.5]);
% 
% %%
% % for i=1:length(Asp)
% %     [maxden maxind]=max(dyn_profile(i).cell_density(54:end,:));
% %     ind=maxind+53;
% %     for j=10:length(ind)
% %         maxind=find(diff(sign(diff(dyn_profile(i).cell_density(:,j))))==-2)+1;
% % 
% %         shift_density(:,j-9)=dyn_profile(i).cell_density(maxind(end)-60:maxind(end)+60,j);
% %     end
% % %     plot(shift_density);hold on;
% %     mean_shift_density(:,i)=mean(shift_density');clear shift_density;
% % end
% % % figure;
% % dx=50/1000;
% % x=0:dx:dx*120;
% % plot(x-50*dx,mean_shift_density*8.4*10^8);
% %%
% % 
% % CWBias_exp=[0,0.0100000000000000,0.0200000000000000,0.0300000000000000,0.0400000000000000,0.0500000000000000,0.0600000000000000,0.0700000000000000,0.0800000000000000,0.0900000000000000,0.100000000000000,0.110000000000000,0.120000000000000,0.130000000000000,0.140000000000000,0.150000000000000,0.160000000000000,0.170000000000000,0.180000000000000,0.190000000000000,0.200000000000000,0.210000000000000,0.220000000000000,0.230000000000000,0.240000000000000,0.250000000000000,0.260000000000000,0.270000000000000,0.280000000000000,0.290000000000000,0.300000000000000,0.310000000000000,0.320000000000000,0.330000000000000,0.340000000000000,0.350000000000000,0.360000000000000,0.370000000000000,0.380000000000000,0.390000000000000,0.400000000000000,0.410000000000000,0.420000000000000,0.430000000000000,0.440000000000000,0.450000000000000,0.460000000000000,0.470000000000000,0.480000000000000,0.490000000000000,0.500000000000000,0.510000000000000,0.520000000000000,0.530000000000000,0.540000000000000,0.550000000000000,0.560000000000000,0.570000000000000,0.580000000000000,0.590000000000000,0.600000000000000,0.610000000000000,0.620000000000000,0.630000000000000,0.640000000000000,0.650000000000000,0.660000000000000,0.670000000000000,0.680000000000000,0.690000000000000,0.700000000000000];
% % P_exp=[0,0.00302823490236943,0.00498640060510877,0.00740654246062159,0.0129576724735464,0.0282364445536801,0.0342746634745129,0.0605832532805051,0.0797351468975508,0.105863629784042,0.144957679208212,0.222281447976504,0.281608928805501,0.373987712933409,0.512065751136908,0.691152325282671,0.919556834551468,1.19485305314468,1.51155357679430,1.88448868628358,2.26163802431774,2.67050119556754,3.04502589236555,3.42156794182481,3.78827243915744,4.11858437055020,4.26770099008759,4.46704480752076,4.50323811436766,4.50544966346908,4.42742856648661,4.31612752213579,4.11316467718262,3.94529362322380,3.65109481682462,3.23848249620281,2.98678093597204,2.69627429599337,2.36626462681502,2.10807472740901,1.95992530735338,1.71658120173055,1.53144380527692,1.43826939988074,1.27455605283792,1.14836389445003,1.01135546697596,0.907225653941910,0.780385038989034,0.685196631521439,0.645362573778435,0.620271119730271,0.579787157033761,0.525204939995039,0.523233085565830,0.466362640857505,0.445131198472857,0.399401650650870,0.399447904195710,0.371516048640824,0.380986255044554,0.393076006538824,0.375160618566594,0.357271202286408,0.341376560968248,0.322917046544370,0.291493170711647,0.311004430859554,0.290374229335908,0.272975724047314,0.162176873861474];
% % dCWBias=0.005;CWBias=0.001:dCWBias:CWBias_exp(end);
% % P=interp1(CWBias_exp,P_exp,CWBias);P=P./sum(P)./dCWBias;
% % figure;
% % plot(CWBias,P); hold on;
% % for i=1:length(Asp)
% % dyn_den_end=final_profile(i).cell_density(:,50+maxind-60:end);
% % P_wave(:,i)=sum(dyn_den_end')./sum(sum(dyn_den_end'))./dCWBias;
% % rel_wave(:,i)= P_wave(:,i)'./P;
% % plot(CWBias,P_wave(:,i));
% % hold on;
% % end
% % %%
% % % % load('C:\Users\XF\Dropbox\paper_wave_diversity\Wave_KSModel_Numerical\test4.mat');
% % % clear all;
% % % Aspind=log(30):(log(40)-log(30))/36:log(40);
% % % % Timeind=log(20):(log(55)-log(20))/100:log(55);
% % % Asp=exp(Aspind);Time=15:6/(length(Asp)-1):25;
% % % tic;
% % % parfor i=1:length(Asp)
% % % %  Asp=[50 100 200];
% % % %  Time=[25 30 35];tic;
% % % %  figure;
% % % %  for i=1:3
% % %        [final_profile(i) dyn_profile]=chemo_oxygen_fun_v2_no_oxg_grad(Asp(i),Time(i)*60);
% % % end
% % %  
% % % save D:\Asp_range_2.mat -v7.3
% % % toc    
% % % figure;
% % % dx=50/1000;
% % % x=0:dx:dx*(length(final_profile(3).density)-1);
% % % for i=1:length(Asp)
% % %     [maxden maxind]=max(final_profile(i).density(50:end));
% % % plot(x-x(maxind+49),final_profile(i).density*8.4*10^8);hold on;
% % % 
% % % end
% % % xlim([-1.5 1.5]);
% %% total cell number
% % load('D:\asp_range.mat');
% % 
% % dx=50;
% % for i=1:length(Asp)
% %      [denmax,denmaxind(i)]=max(final_profile(i).density(200:end));
% %     cell_number(i)=sum(mean_shift_density(:,i))*8.4*10^8*10^-12*600*14*dx;
% % %     plot(final_profile(i).density(200+denmaxind(i)-50:end));
% % %     hold on;
% %     
% % end
% % loglog(Asp,cell_number(1:length(Asp))); hold on;
% % xlim([40 210]);
% % xlabel('Asp (\muM)');ylabel('cell number');
% % %% wave speed
% % % errorbar(asp,wave_speed*60/1000,wave_speed_sd*60/1000,'o');
% % % hold on;
% % for i=1:length(Asp)
% %     peakpos=final_profile(i).peakpos;
% %     time=0:0.5:(length(peakpos)-1)*0.5;
% %     p=polyfit(time(peakpos>7),peakpos(peakpos>7),1);
% %     wavespeed(i)=p(1)*60;
% % end
% % semilogx(Asp, [final_profile.wavespeed]);
% % hold on;
% % semilogx(Asp,  wavespeed)
% % % 
% % %% consumption
% % % exp_con=wave_speed*60/1000.*asp*10^-12*0.6*0.014./total_cells;
% % % plot(total_cells,exp_con,'o');
% % % hold on;
% % sim_con=wavespeed.*Asp*10^-12*0.6*0.014./cell_number;
% % plot(cell_number,sim_con);
% % 
% % % %
% % % 
% % %% fold enrichment
% % CWBias_exp=[0,0.0100000000000000,0.0200000000000000,0.0300000000000000,0.0400000000000000,0.0500000000000000,0.0600000000000000,0.0700000000000000,0.0800000000000000,0.0900000000000000,0.100000000000000,0.110000000000000,0.120000000000000,0.130000000000000,0.140000000000000,0.150000000000000,0.160000000000000,0.170000000000000,0.180000000000000,0.190000000000000,0.200000000000000,0.210000000000000,0.220000000000000,0.230000000000000,0.240000000000000,0.250000000000000,0.260000000000000,0.270000000000000,0.280000000000000,0.290000000000000,0.300000000000000,0.310000000000000,0.320000000000000,0.330000000000000,0.340000000000000,0.350000000000000,0.360000000000000,0.370000000000000,0.380000000000000,0.390000000000000,0.400000000000000,0.410000000000000,0.420000000000000,0.430000000000000,0.440000000000000,0.450000000000000,0.460000000000000,0.470000000000000,0.480000000000000,0.490000000000000,0.500000000000000,0.510000000000000,0.520000000000000,0.530000000000000,0.540000000000000,0.550000000000000,0.560000000000000,0.570000000000000,0.580000000000000,0.590000000000000,0.600000000000000,0.610000000000000,0.620000000000000,0.630000000000000,0.640000000000000,0.650000000000000,0.660000000000000,0.670000000000000,0.680000000000000,0.690000000000000,0.700000000000000];
% % P_exp=[0,0.00302823490236943,0.00498640060510877,0.00740654246062159,0.0129576724735464,0.0282364445536801,0.0342746634745129,0.0605832532805051,0.0797351468975508,0.105863629784042,0.144957679208212,0.222281447976504,0.281608928805501,0.373987712933409,0.512065751136908,0.691152325282671,0.919556834551468,1.19485305314468,1.51155357679430,1.88448868628358,2.26163802431774,2.67050119556754,3.04502589236555,3.42156794182481,3.78827243915744,4.11858437055020,4.26770099008759,4.46704480752076,4.50323811436766,4.50544966346908,4.42742856648661,4.31612752213579,4.11316467718262,3.94529362322380,3.65109481682462,3.23848249620281,2.98678093597204,2.69627429599337,2.36626462681502,2.10807472740901,1.95992530735338,1.71658120173055,1.53144380527692,1.43826939988074,1.27455605283792,1.14836389445003,1.01135546697596,0.907225653941910,0.780385038989034,0.685196631521439,0.645362573778435,0.620271119730271,0.579787157033761,0.525204939995039,0.523233085565830,0.466362640857505,0.445131198472857,0.399401650650870,0.399447904195710,0.371516048640824,0.380986255044554,0.393076006538824,0.375160618566594,0.357271202286408,0.341376560968248,0.322917046544370,0.291493170711647,0.311004430859554,0.290374229335908,0.272975724047314,0.162176873861474];
% % dCWBias=0.002;CWBias=0.01:dCWBias:CWBias_exp(end);
% % P=interp1(CWBias_exp,P_exp,CWBias);P=P./sum(P)./dCWBias;
% % for i=1:length(Asp)
% %     dyn_den_end=final_profile(i).cell_density(:,200+denmaxind(i)-90:end);
% %     P_wave(:,i)=sum(dyn_den_end')./sum(sum(dyn_den_end'))./dCWBias;
% %     rel_wave(:,i)= P_wave(:,i)'./P;
% %  end
% %  surf( CWBias,Asp, rel_wave', 'EdgeColor', 'none');
% %  view(2);
% %  
% %  %%
% %  for i=1:length(Asp)
% %     P_cumsum(:,i)=cumsum(P_wave(:,i))./sum(P_wave(:,i));
% %     [min_cumsum ind_80(i)]=min(abs(P_cumsum(:,i)-0.8));
% %     CW_80(i)=CWBias(ind_80(i));
% %  end
% %  plot(Asp,smooth(CW_80));
% %  hold on;
% %  %%
% %   for i=1:length(Asp)
% %     P_cumsum(:,i)=cumsum(P_wave(:,i))./sum(P_wave(:,i));
% %     [min_cumsum ind_70(i)]=min(abs(P_cumsum(:,i)-0.7));
% %     CW_70(i)=CWBias(ind_70(i));
% %  end
% %  plot(Asp,smooth(CW_70));
% %  hold on;
% %  %%
% %    for i=1:length(Asp)
% %     P_cumsum(:,i)=cumsum(P_wave(:,i))./sum(P_wave(:,i));
% %     [min_cumsum ind_90(i)]=min(abs(P_cumsum(:,i)-0.9));
% %     CW_90(i)=CWBias(ind_90(i));
% %  end
% %  plot(Asp,smooth(CW_90));
% %  
% %  %%
% %  
% %    for i=1:length(Asp)
% %     P_cumsum(:,i)=cumsum(P_wave(:,i))./sum(P_wave(:,i));
% %     [min_cumsum ind_82(i)]=min(abs(P_cumsum(:,i)-0.82));
% %     CW_82(i)=CWBias(ind_82(i));
% %  end
% %  plot(Asp,smooth(CW_82));
%  %% theoritical boundary
% %  Ki=0.3;      %uM
% % Ka=1000;    %uM
% % for i=1:length(Asp)
% %     ASP_profile(:,i)=final_profile(i).asp(:);
% % %     profile1=sum(dyn_profile(i).cell_density(:,:,end-timeind(i)));
% % %     profile2=sum(dyn_profile(i).cell_density(:,:,end-timeind(i)-1));
% % %     [a,ind1]=max(profile1(150:end));
% % %     [a,ind2]=max(profile2(150:end));
% %     instant_wave_speed(i)=final_profile(i).wavespeed/60*1000;
% % end
% %  f=log((1+ASP_profile/Ki)./(1+ASP_profile/Ka));
% %  diff_f=(f(3:end,:)-f(1:end-2,:))/2/dx;
% %  max_diff_f=max(diff_f);
% %  ki_min=smooth(smooth(instant_wave_speed))./max_diff_f';
% %  
% %  for i=1:length(Asp)
% %     options=optimset('TolX',10^-100);
% %     CWBias_max(i)=fzero(@(cwbias)findcwbias(cwbias,ki_min(i)),[0.1,0.8],options);
% %  end
% % z=ones(1,251)*2;
% %  plot3(CWBias_max(Asp>33),Asp(Asp>33),z(Asp>33));
% %  %%
% %  surf( CWBias,smooth([final_profile.wavespeed]), rel_wave', 'EdgeColor', 'none');
% %  view(2);
% %  hold on;
% %  
% %  plot3(CWBias_max(Asp>33),smooth([final_profile(Asp>33).wavespeed]),z(Asp>33));
% % %  %%
% % %  
% %% FigS7B
% %  load('D:\asp_run.mat');
% % x=0:20:20*749;
% % dx=20;
% % for i=1:3
% % cell_density(i,:)=sum(final_profile(i).cell_density);
% % [denmax(i),denmaxind(i)]=max(cell_density(i,300:end));
% % plot((x-(denmaxind(i)-1+300)*dx)/1000,cell_density(i,:)*8.4*10^8);hold on;
% % 
% % end
% % 
% % for i=1:length(Asp)
% %     [maxden maxind]=max(dyn_profile(i).cell_density(100:end,:));
% %     ind=maxind+99;
% %     for j=30:length(ind)
% %         maxind=find(diff(sign(diff(dyn_profile(i).cell_density(:,j))))==-2)+1;
% % 
% %         shift_density(:,j-29)=dyn_profile(i).cell_density(maxind(end)-150:maxind(end)+150,j);
% %     end
% % %     plot(shift_density);hold on;
% %     mean_shift_density(:,i)=mean(shift_density');clear shift_density;
% % end
% % % figure;
% % dx=20/1000;
% % x=0:dx:dx*300;
% % plot(x-150*dx,mean_shift_density*8.4*10^8);
% % 
% % 
% % % %% Fig2insert
% % CWBias_exp=[0,0.0100000000000000,0.0200000000000000,0.0300000000000000,0.0400000000000000,0.0500000000000000,0.0600000000000000,0.0700000000000000,0.0800000000000000,0.0900000000000000,0.100000000000000,0.110000000000000,0.120000000000000,0.130000000000000,0.140000000000000,0.150000000000000,0.160000000000000,0.170000000000000,0.180000000000000,0.190000000000000,0.200000000000000,0.210000000000000,0.220000000000000,0.230000000000000,0.240000000000000,0.250000000000000,0.260000000000000,0.270000000000000,0.280000000000000,0.290000000000000,0.300000000000000,0.310000000000000,0.320000000000000,0.330000000000000,0.340000000000000,0.350000000000000,0.360000000000000,0.370000000000000,0.380000000000000,0.390000000000000,0.400000000000000,0.410000000000000,0.420000000000000,0.430000000000000,0.440000000000000,0.450000000000000,0.460000000000000,0.470000000000000,0.480000000000000,0.490000000000000,0.500000000000000,0.510000000000000,0.520000000000000,0.530000000000000,0.540000000000000,0.550000000000000,0.560000000000000,0.570000000000000,0.580000000000000,0.590000000000000,0.600000000000000,0.610000000000000,0.620000000000000,0.630000000000000,0.640000000000000,0.650000000000000,0.660000000000000,0.670000000000000,0.680000000000000,0.690000000000000,0.700000000000000];
% % P_exp=[0,0.00302823490236943,0.00498640060510877,0.00740654246062159,0.0129576724735464,0.0282364445536801,0.0342746634745129,0.0605832532805051,0.0797351468975508,0.105863629784042,0.144957679208212,0.222281447976504,0.281608928805501,0.373987712933409,0.512065751136908,0.691152325282671,0.919556834551468,1.19485305314468,1.51155357679430,1.88448868628358,2.26163802431774,2.67050119556754,3.04502589236555,3.42156794182481,3.78827243915744,4.11858437055020,4.26770099008759,4.46704480752076,4.50323811436766,4.50544966346908,4.42742856648661,4.31612752213579,4.11316467718262,3.94529362322380,3.65109481682462,3.23848249620281,2.98678093597204,2.69627429599337,2.36626462681502,2.10807472740901,1.95992530735338,1.71658120173055,1.53144380527692,1.43826939988074,1.27455605283792,1.14836389445003,1.01135546697596,0.907225653941910,0.780385038989034,0.685196631521439,0.645362573778435,0.620271119730271,0.579787157033761,0.525204939995039,0.523233085565830,0.466362640857505,0.445131198472857,0.399401650650870,0.399447904195710,0.371516048640824,0.380986255044554,0.393076006538824,0.375160618566594,0.357271202286408,0.341376560968248,0.322917046544370,0.291493170711647,0.311004430859554,0.290374229335908,0.272975724047314,0.162176873861474];
% % dCWBias=0.001;CWBias=0.01:dCWBias:CWBias_exp(end);
% % P=interp1(CWBias_exp,P_exp,CWBias);P=P./sum(P)./dCWBias;
% % figure;
% % 
% % plot(CWBias,P); hold on;
% % for i=1:3
% %  dyn_den_end=final_profile(i).cell_density(:,300+denmaxind(i)-120:end);
% %     P_wave(:,i)=sum(dyn_den_end')./sum(sum(dyn_den_end'))./dCWBias;
% % rel_wave(:,i)= P_wave(:,i)'./P;
% % plot(CWBias,P_wave(:,i));
% % hold on;
% % end
% % figure;
% % plot(CWBias,rel_wave);
% 
% %%
% % for i=1:3
% % plot(CWBias,rel_wave(:,i).*sum(mean_shift_density(:,i))./sum(sum(final_profile(i).cell_density(:,2:end-1))));
% % hold on;
% % end
% % 
% % 
%  %% Fig3D
% dx=20;
% %   final_profile=final_profile(3);
%  cell_density=sum(final_profile.cell_density);
%     [denmax,denmaxind]=max(cell_density(300:end));
%     denmaxind=denmaxind+299;
%     wavedensity=sum(cell_density((200+denmaxind-200):end));
%     cell_number=wavedensity*(dx*6*8.4/10*1.4);
%     
%   dyn_den_end=final_profile.cell_density(:,:);
%   [phenotypemax,phenotypemaxind]=max(dyn_den_end');
%   x=0:dx:(650-1)*dx;
%   x=x-(denmaxind-1)*dx;
%   x=x/1000;
%   yyaxis left;
%   plot(x,sum(dyn_den_end)*8.4*10^8,'DisplayName','cell_density');
%   
%    hold on;
%   CWBias_exp=[0,0.0100000000000000,0.0200000000000000,0.0300000000000000,0.0400000000000000,0.0500000000000000,0.0600000000000000,0.0700000000000000,0.0800000000000000,0.0900000000000000,0.100000000000000,0.110000000000000,0.120000000000000,0.130000000000000,0.140000000000000,0.150000000000000,0.160000000000000,0.170000000000000,0.180000000000000,0.190000000000000,0.200000000000000,0.210000000000000,0.220000000000000,0.230000000000000,0.240000000000000,0.250000000000000,0.260000000000000,0.270000000000000,0.280000000000000,0.290000000000000,0.300000000000000,0.310000000000000,0.320000000000000,0.330000000000000,0.340000000000000,0.350000000000000,0.360000000000000,0.370000000000000,0.380000000000000,0.390000000000000,0.400000000000000,0.410000000000000,0.420000000000000,0.430000000000000,0.440000000000000,0.450000000000000,0.460000000000000,0.470000000000000,0.480000000000000,0.490000000000000,0.500000000000000,0.510000000000000,0.520000000000000,0.530000000000000,0.540000000000000,0.550000000000000,0.560000000000000,0.570000000000000,0.580000000000000,0.590000000000000,0.600000000000000,0.610000000000000,0.620000000000000,0.630000000000000,0.640000000000000,0.650000000000000,0.660000000000000,0.670000000000000,0.680000000000000,0.690000000000000,0.700000000000000];
% P_exp=[0,0.00302823490236943,0.00498640060510877,0.00740654246062159,0.0129576724735464,0.0282364445536801,0.0342746634745129,0.0605832532805051,0.0797351468975508,0.105863629784042,0.144957679208212,0.222281447976504,0.281608928805501,0.373987712933409,0.512065751136908,0.691152325282671,0.919556834551468,1.19485305314468,1.51155357679430,1.88448868628358,2.26163802431774,2.67050119556754,3.04502589236555,3.42156794182481,3.78827243915744,4.11858437055020,4.26770099008759,4.46704480752076,4.50323811436766,4.50544966346908,4.42742856648661,4.31612752213579,4.11316467718262,3.94529362322380,3.65109481682462,3.23848249620281,2.98678093597204,2.69627429599337,2.36626462681502,2.10807472740901,1.95992530735338,1.71658120173055,1.53144380527692,1.43826939988074,1.27455605283792,1.14836389445003,1.01135546697596,0.907225653941910,0.780385038989034,0.685196631521439,0.645362573778435,0.620271119730271,0.579787157033761,0.525204939995039,0.523233085565830,0.466362640857505,0.445131198472857,0.399401650650870,0.399447904195710,0.371516048640824,0.380986255044554,0.393076006538824,0.375160618566594,0.357271202286408,0.341376560968248,0.322917046544370,0.291493170711647,0.311004430859554,0.290374229335908,0.272975724047314,0.162176873861474];
% dCWBias=0.001;CWBias=0.001:dCWBias:CWBias_exp(end);
% P=interp1(CWBias_exp,P_exp,CWBias,'spline');P=P./sum(P)./dCWBias;
%   plot(x,sum(dyn_den_end(CWBias<0.2,:))*8.4*10^8);
%   plot(x,sum(dyn_den_end(abs(CWBias-0.25)<0.05,:))*8.4*10^8);
%   plot(x,sum(dyn_den_end(abs(CWBias-0.35)<0.05,:))*8.4*10^8);
%   plot(x,sum(dyn_den_end(abs(CWBias-0.45)<0.05,:))*8.4*10^8);
%   Ki=2;      %uM
% Ka=1000;    %uM
% N=6;
% ASP_profile=final_profile.asp;
%  f=N*log((1+ASP_profile/Ki)./(1+ASP_profile/Ka));
%  yyaxis right
%  plot(x,f);
%  ylim([0 50]);
%  gradient=((f(3:end)-f(1:end-2))/2/dx)*1000;
%  [maxgrad,maxgradind]=max(gradient);
%  plot((x(1:end-2)+x(3:end))/2,gradient)
% % % 
% %%
% figure;
% yyaxis right;
%  xlim([-3 1]);ylim([0.3 30]);ylabel('gradient mm^{-1}');xlabel('space mm')
% semilogy((x(1:end-2)+x(3:end))/2,gradient);hold on;
% semilogy([-1.56 -1.56],[0.0001 15^10]);ylim([0.3 30]);
%   %%
%  hold on;
% yyaxis left
%    x=0:dx:(650-1)*dx;
% den_profile=final_profile.cell_density(:,1:end);
% [denmax,denmaxind]=max(sum(den_profile(:,200:end)));
% denmaxind=199+denmaxind;
% [phenotypemax,phenotypemaxind]=max(den_profile(:,denmaxind-300:end)');
% x=x-(denmaxind-1)*dx; 
% x=x/1000;
% X_max=x(phenotypemaxind+denmaxind-300);
% CWBias_exp=[0,0.0100000000000000,0.0200000000000000,0.0300000000000000,0.0400000000000000,0.0500000000000000,0.0600000000000000,0.0700000000000000,0.0800000000000000,0.0900000000000000,0.100000000000000,0.110000000000000,0.120000000000000,0.130000000000000,0.140000000000000,0.150000000000000,0.160000000000000,0.170000000000000,0.180000000000000,0.190000000000000,0.200000000000000,0.210000000000000,0.220000000000000,0.230000000000000,0.240000000000000,0.250000000000000,0.260000000000000,0.270000000000000,0.280000000000000,0.290000000000000,0.300000000000000,0.310000000000000,0.320000000000000,0.330000000000000,0.340000000000000,0.350000000000000,0.360000000000000,0.370000000000000,0.380000000000000,0.390000000000000,0.400000000000000,0.410000000000000,0.420000000000000,0.430000000000000,0.440000000000000,0.450000000000000,0.460000000000000,0.470000000000000,0.480000000000000,0.490000000000000,0.500000000000000,0.510000000000000,0.520000000000000,0.530000000000000,0.540000000000000,0.550000000000000,0.560000000000000,0.570000000000000,0.580000000000000,0.590000000000000,0.600000000000000,0.610000000000000,0.620000000000000,0.630000000000000,0.640000000000000,0.650000000000000,0.660000000000000,0.670000000000000,0.680000000000000,0.690000000000000,0.700000000000000];
% P_exp=[0,0.00302823490236943,0.00498640060510877,0.00740654246062159,0.0129576724735464,0.0282364445536801,0.0342746634745129,0.0605832532805051,0.0797351468975508,0.105863629784042,0.144957679208212,0.222281447976504,0.281608928805501,0.373987712933409,0.512065751136908,0.691152325282671,0.919556834551468,1.19485305314468,1.51155357679430,1.88448868628358,2.26163802431774,2.67050119556754,3.04502589236555,3.42156794182481,3.78827243915744,4.11858437055020,4.26770099008759,4.46704480752076,4.50323811436766,4.50544966346908,4.42742856648661,4.31612752213579,4.11316467718262,3.94529362322380,3.65109481682462,3.23848249620281,2.98678093597204,2.69627429599337,2.36626462681502,2.10807472740901,1.95992530735338,1.71658120173055,1.53144380527692,1.43826939988074,1.27455605283792,1.14836389445003,1.01135546697596,0.907225653941910,0.780385038989034,0.685196631521439,0.645362573778435,0.620271119730271,0.579787157033761,0.525204939995039,0.523233085565830,0.466362640857505,0.445131198472857,0.399401650650870,0.399447904195710,0.371516048640824,0.380986255044554,0.393076006538824,0.375160618566594,0.357271202286408,0.341376560968248,0.322917046544370,0.291493170711647,0.311004430859554,0.290374229335908,0.272975724047314,0.162176873861474];
% dCWBias=0.001;CWBias=0.001:dCWBias:CWBias_exp(end);
% P=interp1(CWBias_exp,P_exp,CWBias,'spline');P=P./sum(P)./dCWBias;
% g1=40;g2=40;
% Kd      =  3.06;
% w       =  1.3;
% Yp=g2*Kd./(g2-g1/2+log(1./CWBias-1))-Kd;
% deltaG  =  g1/4-g2/2*(Yp./(Kd+Yp));
% lambdaT   =  w*exp(deltaG);
% lambdaR  =  w*exp(-deltaG);
% CWBias  =  ((lambdaR)./((lambdaT)+(lambdaR)));
% alpha=5.1138;
% a=Yp/alpha;
% Drot=0.062*2;
% v=26*4/pi;%2D run speed 26um/s
% theta=1-0.1564;
% % Diff=v^2*(1-CWBias)./3./(lambdaR*theta+Drot);
% lambdaRdiff=-(lambdaR).*a.*(a-1).*(alpha*g2*Kd)./(2*(Kd+alpha*a).^2); 
% tau=3;%% tau=20.8*exp(-5.86*CWBias); %Diversity in adapation time from Park & Cluzel - Nature (2010)
% N=6;
% %  tau=0;
% ki=theta*v^2*tau.*(lambdaRdiff./(1+lambdaR./lambdaT)./(lambdaR*theta+Drot))./(1+tau.*(theta*lambdaR+Drot))/3;
% ki=ki*10^-6*60;
% semilogy(smooth(smooth(X_max(8:596))),smooth(ki(8:596)));
% ylabel('\chi mm^2/min');
% yyaxis right
% % ASP_profile=dyn_profile(3).Asp(1:end,end);
% %  f=N*log((1+ASP_profile/Ki)./(1+ASP_profile/Ka));
% % 
% % gradient=((f(3:end)-f(1:end-2))/2/dx)*1000;
% % semilogy((x(1:end-2)+x(3:end))/2,gradient);
% % xlim([-3 1]);ylim([0.3 30]);ylabel('gradient mm^{-1}');xlabel('space mm')
% % hold on
% % % 
% % 
% % [maxgrad,maxgradind]=max(gradient);
% % yyaxis right
% % ddx=(x(1:end-2)+x(3:end))/2;
% % plot(ones(1,100)*ddx(maxgradind),0.1:1:99.1);
% yyaxis left
% [cw01 cw01ind]=min(abs(CWBias-0.1));
% plot(X_max(cw01ind),ki(cw01ind),'o');
% [cw02 cw02ind]=min(abs(CWBias-0.2));
% plot(X_max(cw02ind),ki(cw02ind),'o')
% [cw03 cw03ind]=min(abs(CWBias-0.3));
% plot(X_max(cw03ind),ki(cw03ind),'o')
% [cw04 cw04ind]=min(abs(CWBias-0.4));
% plot(X_max(cw04ind),ki(cw04ind),'o');
% yyaxis right
% gradientF=gradient(phenotypemaxind+denmaxind-300-1);
% 
% plot(X_max(cw01ind),gradientF(cw01ind),'o');
% 
% plot(X_max(cw02ind),gradientF(cw02ind),'o')
% 
% plot(X_max(cw03ind),gradientF(cw03ind),'o')
% 
% plot(X_max(cw04ind),gradientF(cw04ind),'o');
% % % 
%   %%
% ind_performance=smooth(smooth(gradientF(13:596).*ki(13:596)));
% 
% 
% plot(smooth(X_max(13:596)),ind_performance);
% hold on;
% dt=0.08;
%  peakpos=final_profile.peakpos;
%  time=0:dt:(length(peakpos)-1)*dt;
%  T=length(peakpos);
% p=polyfit(time(peakpos>7.15),peakpos(peakpos>7.15),1);
% p(1)*60
% x=0:dx:(650-1)*dx;
% x=x-(denmaxind-1)*dx; 
% plot(x/1000,ones(1,length(x)).*p(1)*60,'--');
% xlim([-2 1]);
% 
% plot(X_max(cw01ind),ind_performance(cw01ind),'o');
% 
% plot(X_max(cw02ind),ind_performance(cw02ind),'o')
% 
% plot(X_max(cw03ind),ind_performance(cw03ind),'o')
% 
% plot(X_max(cw04ind),ind_performance(cw04ind),'o');
% ylim([0 0.3]);
% plot(ones(1,100)*ddx(maxgradind),0:1:99);
% %  %% FigS6B
% % dyn_density=sum(dyn_profile(3).cell_density,1);
% % x=0:dx:(750-1)*dx;time=0:dt*60:(length(peakpos)-1)*dt;
% %  surf( x/1000,time/60, squeeze(dyn_density)', 'EdgeColor', 'none');
% %  view(2);
% % xlabel('space mm');
% % ylabel('time min');
% % % 
% % %% FigS6A
% % CWBias=0:0.001:0.99;
% % g1=40;g2=40;
% % Kd      =  3.06;
% % w       =  1.3;
% % Yp=g2*Kd./(g2-g1/2+log(1./CWBias-1))-Kd;
% % deltaG  =  g1/4-g2/2*(Yp./(Kd+Yp));
% % lambdaT   =  w*exp(deltaG);
% % lambdaR  =  w*exp(-deltaG);
% % CWBias  =  ((lambdaR)./((lambdaT)+(lambdaR)));
% % alpha=6;
% % a=Yp/alpha;
% % Drot=0.062*2;
% % v=26*4/pi;%2D run speed 26um/s
% % theta=1-0.1564;
% % Diff=v^2*(1-CWBias)./3./(lambdaR*theta+Drot);
% % lambdaRdiff=-(lambdaR).*a.*(a-1).*(alpha*g2*Kd)./(2*(Kd+alpha*a).^2); 
% % tau=15*exp(-2.2*CWBias); %Diversity in adapation time from Park & Cluzel - Nature (2010)
% % N=6;
% % %  tau=0;
% % ki=theta*v^2*tau.*(lambdaRdiff./(1+lambdaR./lambdaT)./(lambdaR*theta+Drot))./(1+tau.*(theta*lambdaR+Drot))/3;
% % semilogy(CWBias,Diff,'r',CWBias,ki,'g');
% % xlabel('Tumble bias'); ylabel('diff. coef. $$\mu$$, chemo. coef. $$\chi$$ ($$\mu^2$$/s)');
% % 
% % % %% FigS6H
% % % density=squeeze(dyn_density);
% % % for i=1:475
% % %     max_den(i)=density(find(diff(sign(diff(density(:,i))))==-2)+1,i);
% % %     
% % % end
% % % time=dt*60:60*dt:475*dt*60;
% % % plot(time/60,max_den./max_den(125));
% % % 
%  %% simulation with self-attractant 
% % GSfactor=[0 1 10];
% % 
% % betafactor=1;
% % Time=100;
% % Asp=200;
% % for i=1:3
% %     [wavespeed_late(i) final_profile(i) dyn_profile(i) dx]=single_phenotype_selfattractant3(Asp,Time*60,GSfactor(i),betafactor);
% % end
% % dt=2;
% % for i=1:3
% %     density=sum(dyn_profile(i).cell_density,1);
% %     density=squeeze(density);
% %     for j=1:576
% %         ind=find(diff(sign(diff(density(:,j))))==-2)+1;
% %         max_den(j)=density(ind(end),j);
% %     end
% %     time=0:250*dt:575*dt*25;
% %     plot(time/60,smooth(max_den)./max_den(49));hold on;
% % end
% %%
% % Time=120;
% % Asp=200;
% % [final_profile dyn_profile]=chemo_oxygen_fun_v2_no_oxg_grad(Asp,Time*60);
% % density=dyn_profile.cell_density;
% %     density=squeeze(density);
% %     dt=0.5;
% %     for j=1:240
% %         ind=find(diff(sign(diff(density(:,j))))==-2)+1;
% %         max_den(j)=density(ind(end),j);
% %     end
% %     time=dt*60:60*dt:240*dt*60;
% %     plot(time/60,smooth(max_den)./max_den(20));hold on;
% %      %%
% % %  
% %     dx=50;
% %     x=0:dx:799*dx;
% %     cell_density=sum(final_profile(2).cell_density,1);
% %     [denmax,denmaxind]=max(cell_density(50:800));
% %     hold on;
% %       X=x-x(50+denmaxind);
% %      plot(X/1000,cell_density(1:end)*8.4*10^8,'r');
% % %%
% % dx=50;
% % x=0:dx:799*dx;
% % for i=1:1
% %     cell_density=sum(final_profile(i).cell_density,1);
% %     
% %     [denmax,denmaxind(i)]=max(cell_density(50:800));
% % %     cell_number(i)=sum(final_profile(i).density(200+denmaxind(i)-50:end))*10^9*10^-12*600*14*dx;
% %     X=x-x(50+denmaxind(i));
% %      plot(X/1000,cell_density(1:end));
% %      hold on;
% %     
% % end
% 
% 
% %%
% 
% % % % % % % 
% % % % % %%
% % %%%%%%%Figure 4A
% % Ki=0.3;      %uM
% % Ka=1000;    %uM
% % den_profile=dyn_profile.cell_density(:,1:end,end-100);
% % ASP_profile=dyn_profile.Asp(1:end,end-100);
% % [denmax(i),denmaxind(i)]=max(sum(den_profile(:,:)));
% % x=0:dx:(750-1)*dx;
% %   x=x-(denmaxind-1)*dx;
% % x=x/1000;
% % plot(x,sum(den_profile)*8.5*10^8,'DisplayName','cell_density');
% % % figure;subplot(2,1,1);
% % plot(x,ASP_profile);
% % f=log((1+ASP_profile/Ki)./(1+ASP_profile/Ka));
% % %  
% % % subplot(2,1,2);
% % gradientf=(-f(1:end-2)+f(3:end))/2/dx;
% % plot((x(1:end-2)+x(3:end))/2,smooth(gradientf)*1000);
% % [gradientfmax,gradientfmaxind]=max(gradientf*1000);
% % % hold on;
% % % plot(x(gradientfmaxind),gradientfmax,'ro');
% % % plot(x(gradientfmaxind)*ones(length(gradientf),1),gradientf,'r--');
% % % %  %%
% % %   [phenotypemax,phenotypemaxind]=max(den_profile(:,denmaxind-300:end)');
% % %  X_max=x(phenotypemaxind+denmaxind-300);
% % % plot(smooth(X_max(10:542)),smooth(smooth(ki(10:542))));
% % % % hold on;
% % % % plot(ones(1,length(ki))*x(gradientfmaxind+den_profile_maxind-200),ki,'r--');
% % % % hold on;
% % % % plot(x,ki(503)*ones(1,length(x)),'r--');
% % % % figure;
% % % % plot(CWBias(13:503),X_max(13:503))
% % % % %%
% % % gradientf=(-f(1:end-2)+f(3:end))/2/dx;
% % % gradient_F=gradientf(phenotypemaxind+denmaxind-300-1);
% % % plot(smooth(X_max(10:542)),smooth(smooth(gradient_F(10:542)'.*ki(10:542))));
% % % % [a, ind1]=max(profile1);
% % % % [a,ind2]=max(profile2);
% % % % 
% % % %  hold on;
% % % % 
% % % % plot(X_max(13:503),ones(1,length(X_max(13:503))).*(ind1-ind2)*dx/40,'--');
% % % % ylim([0 6]);
% % % % figure;
% % % % plot(smooth(X_max(13:503)),smooth(gradient_F(13:503)));
% % % % 
% % % % %% x vs mean ki
% % % % for i=1:length(x)
% % % %     meanki(:,i)=sum(ki'.*den_profile(:,i)./sum(den_profile(:,:)')');
% % % % end
% % % % 
% % % % %%
% % % %  %%
% % % % % %%%%%Figure 4B
% % % % % [phenotypemax,phenotypemaxind]=max(den_profile(:,den_profile_maxind-200:den_profile_maxind+100)');
% % % % % plot(CWBias,(x(phenotypemaxind+den_profile_maxind-200)));
% % % % % hold on;
% % % % % plot(CWBias,ones(1,length(CWBias))*x(gradientfmaxind+den_profile_maxind-100),'r--');
% % % % % hold on;
% % % % % plot(CWBias(gradientfmaxind+den_profile_maxind)*ones(1,length(x)),x,'r--');
% % % % % for i=1:length(CWBias)
% % % % %    phenotypemeanpos(i)=sum(den_profile(i,den_profile_maxind-200:end).*x(den_profile_maxind-200:end))./sum(den_profile(i,den_profile_maxind-200:end));
% % % % % end
% % % % % 
% % % % % plot(phenotypemeanpos,CWBias);
% % % % %  %%
% % % % % %%%%%%Figure 4C
% % % % % [gradmax,gradmaxind]=max(gradientf(phenotypemaxind));
% % % % % figure;
% % % % % subplot(2,1,1);
% % % % % hold on;
% % % % % plot(CWBias(1:gradmaxind),smooth(gradientf(phenotypemaxind(1:gradmaxind)))');
% % % % % subplot(2,1,2);
% % % % % plot(CWBias(1:gradmaxind),ki(1:gradmaxind));
% % % % % figure;
% % % % % plot(x(den_profile_maxind-100:den_profile_maxind+100)',sum(den_profile(:,den_profile_maxind-100:den_profile_maxind+100)));
% % % % % % %% Figure 4D
% % % % % figure;
% % % % % 
% % % profile1=sum(dyn_profile(end).cell_density(:,:,end));
% % % profile2=sum(dyn_profile(end).cell_density(:,:,end-200));
% % % [a, ind1]=max(profile1);
% % % [a,ind2]=max(profile2);
% % % % % plot(CWBias(1:gradmaxind),smooth(smooth(gradientf(phenotypemaxind((1:gradmaxind))))'.*ki((1:gradmaxind))));
% % % % %  hold on;
% % % % % 
% % % plot(x,ones(1,length(x)).*(ind1-ind2)*dx/800,'--');
% % % % ylim([0 6]);
%  %% 
% %% 
% 
% %  g1=40;g2=40;
% % Kd      =  3.06;
% % w       =  1.3;
% % Yp=g2*Kd./(g2-g1/2+log(1./CWBias-1))-Kd;
% % deltaG  =  g1/4-g2/2*(Yp./(Kd+Yp));
% % lambdaT   =  w*exp(deltaG);
% % lambdaR  =  w*exp(-deltaG);
% % CWBias  =  ((lambdaR)./((lambdaT)+(lambdaR)));
% % alpha=6;
% % a=Yp/alpha;
% % Drot=0.062*2;
% % tauR=1./(Drot+lambdaR);
% % v=26*4/pi;%2D run speed 26um/s
% % theta=1-0.1564;
% % Diff=v^2*(1-CWBias).*tauR/3/theta;
% % lambdaRdiff=-(lambdaR).*a.*(a-1).*(alpha*g2*Kd)./(2*(Kd+alpha*a).^2); 
% % tau=20.8*exp(-5.86*CWBias); 
% % N=6;
% % Ki=0.3;      %uM
% % Ka=1000;    %uM
% % 
% % %  tau=0;
% % ki=theta*v^2*tau.*(lambdaRdiff./(1+lambdaR./lambdaT)./(lambdaR*theta+Drot))./(1+tau.*(theta*lambdaR+Drot))/3;
% % cell_density=sum(final_profile.cell_density);
% %     [denmax,denmaxind]=max(cell_density(300:end));
% %     denmaxind=denmaxind+299;
% %     cell_number=wavedensity*(dx*6*8.4/10*1.4);
% %     
% %   dyn_den_end=dyn_profile(3).cell_density(:,denmaxind-300:end,end);
% %   [phenotypemax,phenotypemaxind]=max(dyn_den_end');
% %   x=0:dx:(650-1)*dx;
% %   
% %  ASP_profile=dyn_profile(3).Asp(denmaxind-300:end,end);
% %   x=0:dx:(length(ASP_profile)-1)*dx;
% %  f=N*log((1+ASP_profile/Ki)./(1+ASP_profile/Ka));
% % 
% % gradient=((f(3:end)-f(1:end-2))/2/dx);
% % plot(CWBias,smooth(gradient(phenotypemaxind)'.*ki));
% % hold on;
% %  dt=0.08;
% %  peakpos=final_profile.peakpos;
% %  time=0:dt:(length(peakpos)-1)*dt;
% %  T=length(peakpos);
% % p=polyfit(time(end-round(T/100):end),peakpos(end-round(T/100):end),1);
% % p(1)*1000
% % x=0:dx:(650-1)*dx;
% % x=x-(denmaxind-1)*dx; 
% % plot(x/1000,ones(1,length(x)).*p(1)*1000,'--');
%  
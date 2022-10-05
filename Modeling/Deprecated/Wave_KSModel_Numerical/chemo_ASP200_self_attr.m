%% simulation with self-attractant 
GSfactor=[0 0.5];%GSfactor=0.18;

betafactor=0.2;
Time=60;
Asp=200;
for i=1:length(GSfactor)
    [final_profile(i) wavespeed(i) dyn_profile(i) dx]=chemo_self_attractant2(Asp,Time*60,GSfactor(i),betafactor);
end

%%
% GSfactor=[0 0.25 1];%GSfactor=0.18;
% 
% betafactor=0.025;
% Time=30;
% Asp=200;
% for i=1:length(GSfactor)
%     [wavespeed_late(i) final_profile(i) dyn_profile(i) dx]=single_phenotype_selfattractant4(Asp,Time*60,GSfactor(i),betafactor);
% end
% dt=0.5;
%% density profile
dt=0.08;
  co=jet(length(GSfactor));
  dx=20/1000;
   
% for i=1:length(GSfactor)
%     cell_density=sum(final_profile(i).cell_density);
%      x=0:dx:dx*(length(cell_density)-1);
%      [maxden maxind]=max(cell_density(48:end));
%      X_pos(:,i)=x-x(maxind+47);
%      Density_shift(:,i)=cell_density(1:end)*8.4*10^8;
%     
%      plot(x-x(maxind+47),cell_density(1:end)*8.4*10^8,'Color',co(i,:));hold on;xlim([-1.5 1.5]);
% end

% dt=0.5;
%   co=jet(length(GSfactor));
%   dx=50/1000;
%    
% for i=1:length(GSfactor)
%     cell_density=sum(final_profile(i).cell_density);
%      x=0:dx:dx*(length(cell_density)-1);
%      [maxden maxind]=max(cell_density(48:end));
%      X_pos(:,i)=x-x(maxind+47);
%      Density_shift(:,i)=cell_density(1:end)*8.4*10^8;
%      plot(x-x(maxind+47),cell_density(1:end)*8.4*10^8,'Color',co(i,:));hold on;xlim([-1.5 1.5]);
% end
%%

for i=1:length(GSfactor)
    density=sum(dyn_profile(i).cell_density,1);
    density=squeeze(density);
    cell_density=sum(final_profile(i).cell_density);
    x=0:dx:dx*(length(cell_density)-1);
    for j=1:120
         [maxden(j) maxind]=max(density(120:end,j));
        max_den(j,i)=density(maxind+119,j);
        x_pos=x-x(maxind+119);
         Cell_Number(j,i)=sum(density(x_pos>-1.5&x_pos<1.5,j));
    end
    max_den_norm(:,i)=max_den(:,i)./max_den(10,i);
     Cell_Number_norm(:,i)=Cell_Number(:,i)./Cell_Number(10,i);
%     time=dt*25:25*dt:576*dt*25;
%     plot(time/60,smooth(max_den)./max_den(48));hold on;
end
figure;
plot([1:120]/2,max_den_norm,'DisplayName','max_den')
% 
% Asp=200;Time=120;
% i=1;
% TB=0.3357;
% [wavespeed_late(i) final_profile(i) dyn_profile(i) dx]=single_phenotype_simple(Asp(i),Time(i)*60,TB);
% %%
% Cell_density=dyn_profile.cell_density(:,:,:);
% for i=1:36
% cellden=sum(Cell_density(:,:,i));
% peaks=find(diff(sign(diff(cellden)))==-2)+1;peak0(i)=peaks(end);
% peakden0(i)=cellden(peak0(i));
% end   
% time=0:8:(900-1)*8;
% plot(time,peakden0)
% % %%
% % clear all
% Asp=200;Time=120;
% i=1;
% TB=0.3357;
% [wavespeed_late(i) final_profile(i) dyn_profile(i) dx]=chemo_single_phenotype(Asp(i),Time(i)*60,TB);
% %%
% Cell_density=dyn_profile.cell_density(:,:,:);
% for i=1:900
% cellden=sum(Cell_density(:,:,i));
% peak2(i)=find(diff(sign(diff(cellden)))==-2)+1;
% peakden2(i)=cellden(peak2(i));
% end   
% time=0:8:(900-1)*8;
% hold on;
% plot(time,peakden2)
% T=time/60;
% peakden10=peakden1(T>10);
% peakden20=peakden2(T>10);
% peakden10=peakden10(1);
% peakden20=peakden20(1);
% figure;
% plot(T(T>10),peakden1(T>10)/peakden10,T(T>10),peakden2(T>10)/peakden20)
%% self-attractant
%%
Asp=200;Time=120;
i=1;
TB=0.3357;
GSfacter=10^-2;
betafacter=1;
% figure;
[wavespeed_late(i) final_profile(i) dyn_profile(i) dx]=single_phenotype_selfattractant(Asp(i),Time(i)*60,TB,GSfacter,betafacter);
Cell_density=dyn_profile.cell_density(:,:,:);
for i=1:900
cellden=sum(Cell_density(:,:,i));
peaks=find(diff(sign(diff(cellden)))==-2)+1;peak1(i)=peaks(end);
peakden1(i)=cellden(peak1(i));
end   
time=0:8:(900-1)*8;
% figure;
plot(time,peakden1)
%%
Asp=200;Time=120;
i=1;
TB=0.3357;
GSfacter=10^-1;
betafacter=1;

[wavespeed_late(i) final_profile(i) dyn_profile(i) dx]=single_phenotype_selfattractant(Asp(i),Time(i)*60,TB,GSfacter,betafacter);
Cell_density=dyn_profile.cell_density(:,:,:);
for i=1:900
cellden=sum(Cell_density(:,:,i));
peaks=find(diff(sign(diff(cellden)))==-2)+1;peak2(i)=peaks(end);
peakden2(i)=cellden(peak2(i));
end     
time=0:8:(900-1)*8;
% figure;
hold on
plot(time,peakden2);

%%
Asp=200;Time=120;
i=1;
TB=0.3357;
GSfacter=1;
betafacter=1;

[wavespeed_late(i) final_profile(i) dyn_profile(i) dx]=single_phenotype_selfattractant(Asp(i),Time(i)*60,TB,GSfacter,betafacter);
Cell_density=dyn_profile.cell_density(:,:,:);
for i=1:900
cellden=sum(Cell_density(:,:,i));
peaks=find(diff(sign(diff(cellden)))==-2)+1;peak3(i)=peaks(end);
peakden3(i)=cellden(peak3(i));
end   
time=0:8:(900-1)*8;

plot(time,peakden3)
%%
Asp=200;Time=120;
i=1;
TB=0.3357;
GSfacter=10;
betafacter=1;
figure;
[wavespeed_late(i) final_profile(i) dyn_profile(i) dx]=single_phenotype_selfattractant(Asp(i),Time(i)*60,TB,GSfacter,betafacter);
Cell_density=dyn_profile.cell_density(:,:,:);
for i=1:900
cellden=sum(Cell_density(:,:,i));
peaks=find(diff(sign(diff(cellden)))==-2)+1;
peak4(i)=peaks(end);
peakden4(i)=cellden(peak4(i));
end   
time=0:8:(900-1)*8;

plot(time,peakden4);
%%
Asp=200;Time=120;
i=1;
TB=0.3357;
GSfacter=100;
betafacter=1;
figure;
[wavespeed_late(i) final_profile(i) dyn_profile(i) dx]=single_phenotype_selfattractant(Asp(i),Time(i)*60,TB,GSfacter,betafacter);
Cell_density=dyn_profile.cell_density(:,:,:);
for i=1:900
cellden=sum(Cell_density(:,:,i));
peaks=find(diff(sign(diff(cellden)))==-2)+1;
peak5(i)=peaks(end);
peakden5(i)=cellden(peak5(i));
end   
time=0:8:(900-1)*8;

plot(time,peakden4);
%%
Asp=200;Time=120;
i=1;
TB=0.3357;
GSfacter=1000;
betafacter=1;
figure;
[wavespeed_late(i) final_profile(i) dyn_profile(i) dx]=single_phenotype_selfattractant(Asp(i),Time(i)*60,TB,GSfacter,betafacter);
Cell_density=dyn_profile.cell_density(:,:,:);
for i=1:900
cellden=sum(Cell_density(:,:,i));
peaks=find(diff(sign(diff(cellden)))==-2)+1;
peak6(i)=peaks(end);
peakden6(i)=cellden(peak6(i));
end   
time=0:8:(900-1)*8;

plot(time,peakden4);
%%
T=time/60; 
 peakden10=smooth(peakden1(T>10));
 peakden20=smooth(peakden2(T>10));
 peakden30=smooth(peakden3(T>10));
 peakden40=smooth(peakden4(T>10));
 peakden50=smooth(peakden5(T>10));
 peakden60=smooth(peakden6(T>10));
peakden10=peakden10(1);
peakden20=peakden20(1);
peakden30=peakden30(1);
peakden40=peakden40(1);
peakden50=peakden50(1);
peakden60=peakden60(1);
figure;
plot(T(T>10),peakden1(T>10)/peakden10,T(T>10),peakden2(T>10)/peakden20,T(T>10),peakden3(T>10)/peakden30,T(T>10),peakden4(T>10)/peakden40,T(T>10),peakden5(T>10)/peakden50,T(T>10),peakden6(T>10)/peakden60);
% 
% 
%%
X=0:dx:1999*dx;
plot((X-X(peak4(end)))/1000,sum(final_profile.cell_density(:,:))*0.84);xlim([-1.8 1.3]);ylim([0,8]);

%%
Asp=200;Time=120;
i=1;
TB=0.3357;
GSfacter=0;
betafacter=1;
% figure;
[wavespeed_late(i) final_profile(i) dyn_profile(i) dx]=single_phenotype_selfattractant(Asp(i),Time(i)*60,TB,GSfacter,betafacter);
Cell_density=dyn_profile.cell_density(:,:,:);
for i=1:900
cellden(i,:)=sum(Cell_density(:,:,i));
peaks=find(diff(sign(diff(cellden)))==-2)+1;peak1(i)=peaks(end);
peakden1(i)=cellden(peak1(i));
end   
time=0:8:(900-1)*8;
% % figure;
% % plot(time,peakden1)
X=0:dx:1999*dx;
surf(time,X,cellden*0.84);
% plot((X-X(peak1(end)))/1000,sum(final_profile.cell_density(:,:))*0.84);xlim([-1.8 1.3]);ylim([0,8]);
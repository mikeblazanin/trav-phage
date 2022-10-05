
Asp=200;Time=120;
i=1;
TB=0.3357;
[wavespeed_late(i) final_profile(i) dyn_profile(i) dx]=single_phenotype_simple(Asp(i),Time(i)*60,TB);
%%
Cell_density=dyn_profile.cell_density(:,:,:);
for i=1:36
cellden=sum(Cell_density(:,:,i));
peaks=find(diff(sign(diff(cellden)))==-2)+1;peak0(i)=peaks(end);
peakden0(i)=cellden(peak0(i));
end   
time=0:8:(900-1)*8;
plot(time,peakden0)

 
 %%
 Asp=200;Time=120;
i=1;
TB=0.3357;
[wavespeed_late(i) final_profile(i) dyn_profile(i) dx]=chemo_single_phenotype(Asp(i),Time(i)*60,TB);

%%
Cell_density=dyn_profile.cell_density(:,:,:);
for i=1:900
cellden=sum(Cell_density(:,:,i));
peak2(i)=find(diff(sign(diff(cellden)))==-2)+1;
peakden2(i)=cellden(peak2(i));
end   
X=0:dx:1999*dx; 
plot((X-X(peak2(end)))/1000,sum(final_profile.cell_density(:,:))*0.84);xlim([-1.8 1.3]);ylim([0,8]);
time=0:8:(900-1)*8;
hold on;
plot(time,peakden2)
T=time/60;
peakden10=peakden1(T>10);
peakden20=peakden2(T>10);
peakden10=peakden10(1);
peakden20=peakden20(1);
figure;
plot(T(T>10),peakden1(T>10)/peakden10,T(T>10),peakden2(T>10)/peakden20)
%%
% clear all
Asp=200;Time=120;
i=1;
TB=0.3357;
[wavespeed_late(i) final_profile(i) dyn_profile(i) dx]=chemo_single_phenotype(Asp(i),Time(i)*60,TB);
%%
Cell_density=dyn_profile.cell_density(:,:,:);
for i=1:900
cellden=sum(Cell_density(:,:,i));
peak2(i)=find(diff(sign(diff(cellden)))==-2)+1;
peakden2(i)=cellden(peak2(i));
end   
time=0:8:(900-1)*8;
hold on;
plot(time,peakden2)
T=time/60;
peakden10=peakden1(T>10);
peakden20=peakden2(T>10);
peakden10=peakden10(1);
peakden20=peakden20(1);
figure;
plot(T(T>10),peakden1(T>10)/peakden10,T(T>10),peakden2(T>10)/peakden20)



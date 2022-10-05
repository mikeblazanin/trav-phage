load('D:\Dropbox\Yale\paper_wave_diversity\manuscript\Figures\final_figure_for_revision\Fig4\fig4f_stat.mat');

parfor i=1:length(data)
    P_exp1=data(i).TB_pdf1;
    P_exp2=data(i).TB_pdf2;
    CWBias=0:0.0101:1;
    [final_profile wavespeed dyn_profile]=chemo_two_population(CWBias,P_exp1,P_exp2);
    result(i).final=final_profile;
    result(i).wavespeed=wavespeed;
    result(i).dyn_profile=dyn_profile;
  
    
end
%%
load('D:\Dropbox\Yale\paper_wave_diversity\manuscript\Figures\final_figure_for_revision\Fig4\fig4f_stat.mat');

for i=1:length(data)
   cell_den=result(i).final.cell_density;
   cell_den_total=sum(cell_den);
   [den_max_total,den_max_total_ind]=max(cell_den_total(51:end));
   den_max_total_ind=den_max_total_ind+50;
   cell_den1=sum(cell_den(1:198,:));
   cell_den2=sum(cell_den(199:end,:));
   [den_max1,den_max_ind1]=max(cell_den1(den_max_total_ind-50:den_max_total_ind+50));
   [den_max2,den_max_ind2]=max(cell_den2(den_max_total_ind-50:den_max_total_ind+50));
   peak_dist(i)=(den_max_ind1-den_max_ind2)*50/1000;
   tb_diff(i)=-data(i).meanTB1+data(i).meanTB2;
   
end
plot(tb_diff,peak_dist,'o')
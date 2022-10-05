load('D:\Dropbox\Yale\paper_wave_diversity\others\Adam_sims_phenotypes\mu_chi_tb.mat');
for i=1:40
    tbi=0+(i-1)*0.025;tbf=tbi+0.025*i;
    chi_mean(i)=mean(chi((tb>tbi&(tb<tbf))));
    mu_mean(i)=mean(mu((tb>tbi&(tb<tbf)&~isnan(mu))));
end
tb_bin=0.0125:0.025:1;
plot(tb_bin,chi_mean,tb_bin,mu_mean)
plot(tb_bin,chi_mean./mu_mean,'o');
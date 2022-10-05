load('C:\Users\jl2345\Dropbox (emonetlab)\users\junjiajia_long\analysis\wave\Asp_Oxy_HillConc_BellTB_NLrec\chemo_oxygen_fun_v2_no_oxg_grad_200_7200_30_30000_0.5_50\chemo_oxygen_fun_v2_no_oxg_grad_200_7200_30_30000_0.5_50_analysis.mat')
f = P_wave(:,100)';
save('custom_TB_distribution', 'f', 'CWBias')
load('wavespeed_exp_data.mat');
figure;
plot(wavespeed_exp_data.cell_number(wavespeed_exp_data.asp==0.05),wavespeed_exp_data.consumption(wavespeed_exp_data.asp==0.05),'ro');
hold on;
plot(wavespeed_exp_data.cell_number(wavespeed_exp_data.asp==0.1),wavespeed_exp_data.consumption(wavespeed_exp_data.asp==0.1),'bo');
plot(wavespeed_exp_data.cell_number(wavespeed_exp_data.asp==0.2),wavespeed_exp_data.consumption(wavespeed_exp_data.asp==0.2),'go');
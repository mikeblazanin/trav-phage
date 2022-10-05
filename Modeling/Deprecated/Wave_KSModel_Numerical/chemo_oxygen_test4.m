% load('C:\Users\XF\Dropbox\paper_wave_diversity\Wave_KSModel_Numerical\test4.mat');
% Asp=40:5:230;Time=20:35/(length(Asp)-1):55;tic;
% for i=1:length(Asp)
%       [wavespeed_late(i) final_profile(i) dyn_profile(i) dx]=chemo_oxygen_fun_v2_no_oxg_grad(Asp(i),Time(i)*60);
% 
% end
clear all

Asp=200;
Time=40;tic;
i=1;
[wavespeed_late(i) final_profile(i) dyn_profile(i) dx]=chemo_oxygen_fun_cheZ(Asp(i),Time(i)*60);
save atc4_atc6.mat

clear all
figure;
Asp=200;
Time=40;tic;
i=1;
[wavespeed_late(i) final_profile(i) dyn_profile(i) dx]=chemo_oxygen_fun_cheZ2(Asp(i),Time(i)*60);
save atc2_atc40.mat
%%
% for i=1:length(Asp)
%     
%       cell_density(i,:)=sum(final_profile(i).cell_density);
%     [denmax(i),denmaxind(i)]=max(cell_density(i,150:end));
%    
%     timeind(i)=0;
%     if denmaxind(i)>60
%         timeind(i)=1;
%         cell_density(i,:)=sum(dyn_profile(i).cell_density(:,1:end,end-timeind(i)));
%         [denmax(i),denmaxind(i)]=max(cell_density(i,150:end));
%     end
%      if denmaxind(i)>60
%         timeind(i)=2;
%         cell_density(i,:)=sum(dyn_profile(i).cell_density(:,1:end,end-timeind(i)));
%         [denmax(i),denmaxind(i)]=max(cell_density(i,150:end));
%      end
%      wavedensity(i)=sum(cell_density(i,(150+denmaxind(i)-110):end));
%     cell_number(i)=wavedensity(i)*(dx*6*8.4/10*1.4);
%    
%     plot(cell_density(i,(150+denmaxind(i)-110):end));
% hold on;
% end
% figure;
% plot(Asp(Asp>40&Asp<210),(cell_number(Asp>40&Asp<210)));
% figure;
% plot(cell_number(Asp>40&Asp<210),Asp(Asp>40&Asp<210).*wavespeed_late(Asp>40&Asp<210)./cell_number(Asp>40&Asp<210)*10^-18*600*14*60/1000,'r-');
% 
% 
% %  plot(Asp(Asp>40&Asp<210),cell_number(Asp>40&Asp<210));
% %  xlim([20,250]);
%  %%
%   dCWBias=0.0025; CWBias=0.01:dCWBias:0.99;P = lognpdf(CWBias,log(0.295),0.308)/sum(lognpdf(CWBias,log(0.295),0.308))/dCWBias;
%  for i=1:length(Asp)
%     dyn_den_end=dyn_profile(i).cell_density(:,150+denmaxind(i)-60:150+denmaxind(i)+60,end-timeind(i));
%     P_wave(:,i)=sum(dyn_den_end')./sum(sum(dyn_den_end'))./dCWBias;
%     rel_wave(:,i)= P_wave(:,i)'./P;
%  end
%  surf( Asp,CWBias, rel_wave, 'EdgeColor', 'none');
%  view(2);
 
%  
% %%
% Ki=0.3;      %uM
% Ka=1000;    %uM
% 
% 
% %  O_profile=dyn_profile(i).O((200+denmaxind(i)-100):(200+denmaxind(i)+100),end);
% 
% for i=1:length(Asp)
%     ASP_profile(:,i)=dyn_profile(i).Asp(:,end-timeind(i));
%     profile1=sum(dyn_profile(i).cell_density(:,:,end-timeind(i)));
%     profile2=sum(dyn_profile(i).cell_density(:,:,end-timeind(i)-1));
%     [a,ind1]=max(profile1(150:end));
%     [a,ind2]=max(profile2(150:end));
%     instant_wave_speed(i)=(ind1-ind2)*dx/120;
% end
%  f=log((1+ASP_profile/Ki)./(1+ASP_profile/Ka));
%  diff_f=(f(3:end,:)-f(1:end-2,:))/2/dx;
%  max_diff_f=max(diff_f);
%  ki_min=smooth(smooth(instant_wave_speed))./max_diff_f';
%  
%  parfor i=1:length(Asp)
%     options=optimset('TolX',10^-100);
%     CWBias_max(i)=fzero(@(cwbias)findcwbias(cwbias,ki_min(i)),[0.1,0.7],options);
% end
 %%
%  Asp=30:5:300;Time=20:35/(length(Asp)-1):55;tic;
%  i=35;
% %  [wavespeed(i) final_profile(i) dyn_profile(i) dx]=chemo_oxygen_fun_v2_finebin_no_oxg_grad(Asp(i),Time(i)*60);
%  cell_density(i,:)=sum(final_profile(i).cell_density);
%     [denmax(i),denmaxind(i)]=max(cell_density(i,200:end));
%     wavedensity(i)=sum(cell_density(i,(200+denmaxind(i)-100):end));
%     cell_number(i)=wavedensity(i)*(dx*6*8.4/10*1.4);
% %   
%    %%
% %  
%   dyn_den_end=dyn_profile(i).cell_density(:,(200+denmaxind(i)-150):(200+denmaxind(i)+150),end);
%   [phenotypemax,phenotypemaxind]=max(dyn_den_end');
%   x=0:dx:(length((200+denmaxind(i)-150):(200+denmaxind(i)+150))-1)*dx;
%   plot(x,sum(dyn_den_end)*8.5*10^8./sum(sum(dyn_den_end)*8.5*10^8),'DisplayName','cell_density')
%   hold on;
%    dCWBias=0.001; CWBias=0.01:dCWBias:0.99;P = lognpdf(CWBias,log(0.295),0.308)/sum(lognpdf(CWBias,log(0.295),0.308))/dCWBias;
%   plot(x,dyn_den_end(abs(CWBias-0.1)<0.0009,:)*8.5*10^8/sum(dyn_den_end(abs(CWBias-0.1)<0.0009,:)*8.5*10^8));
%   plot(x,dyn_den_end(abs(CWBias-0.2)<0.0009,:)*8.5*10^8/sum(dyn_den_end(abs(CWBias-0.2)<0.0009,:)*8.5*10^8));
%   plot(x,dyn_den_end(abs(CWBias-0.3)<0.0009,:)*8.5*10^8/sum(dyn_den_end(abs(CWBias-0.3)<0.0009,:)*8.5*10^8));
%   plot(x,dyn_den_end(abs(CWBias-0.4)<0.0009,:)*8.5*10^8/sum(dyn_den_end(abs(CWBias-0.4)<0.0009,:)*8.5*10^8));
% %%
% figure;
%   hold on;
% %   dCWBias=0.001; CWBias=0.01:dCWBias:0.99;P = lognpdf(CWBias,log(0.295),0.308)/sum(lognpdf(CWBias,log(0.295),0.308))/dCWBias;
%     plot(x,sum(dyn_den_end)*8.5*10^8,'DisplayName','cell_density')
%   plot(x,sum(dyn_den_end(abs(CWBias-0.15)<0.05,:))*8.5*10^8);
%   plot(x,sum(dyn_den_end(abs(CWBias-0.25)<0.05,:))*8.5*10^8);
%   plot(x,sum(dyn_den_end(abs(CWBias-0.35)<0.05,:))*8.5*10^8);
%   plot(x,sum(dyn_den_end(abs(CWBias-0.45)<0.05,:))*8.5*10^8);
% %  plot(x,(sum(dyn_den_end(abs(CWBias-0.15)<0.05,:))+sum(dyn_den_end(abs(CWBias-0.25)<0.05,:))+sum(dyn_den_end(abs(CWBias-0.35)<0.05,:)))*8.5*10^8);
%  %%
% %   
% Ki=0.3;      %uM
% Ka=1000;    %uM
% KiO1=40;
% KiO2=150;
% KaO1=330;
% KaO2=10000;
%  ASP_profile=dyn_profile(i).Asp((200+denmaxind(i)-150):(200+denmaxind(i)+150),end);
% %  O_profile=dyn_profile(i).O((200+denmaxind(i)-100):(200+denmaxind(i)+100),end);
%  f=log((1+ASP_profile/Ki)./(1+ASP_profile/Ka));
% %  fO=(log((1+O_profile/KiO1)./(1+O_profile/KaO1))-log((1+O_profile/KiO2)./(1+O_profile/KaO2)))/5000;
% %  plot((x(1:end-2)+x(3:end))/2,((f(3:end)-f(1:end-2))/2/dx+(fO(3:end)-fO(1:end-2))/2/dx)./max((f(3:end)-f(1:end-2))/2/dx+(fO(3:end)-fO(1:end-2))/2/dx)*5*10^9);
% plot(x,ASP_profile./200*8*10^9);
%  g1=40;g2=40;
% Kd      =  3.06;
% w       =  1.3;
% Yp=g2*Kd./(g2-g1/2+log(1./CWBias-1))-Kd;
% deltaG  =  g1/4-g2/2*(Yp./(Kd+Yp));
% lambdaT   =  w*exp(deltaG);
% lambdaR  =  w*exp(-deltaG);
% CWBias  =  ((lambdaR)./((lambdaT)+(lambdaR)));
% alpha=6;
% a=Yp/alpha;
% Drot=0.062*2;
% tauR=1./(Drot+lambdaR);
% v=25*4/pi;%2D run speed 26um/s
% theta=1-0.33;
% Diff=v^2*(1-CWBias).*tauR/3/theta;
% lambdaRdiff=-(lambdaR).*a.*(a-1).*(alpha*g2*Kd)./(2*(Kd+alpha*a).^2); 
% tau=20*exp(-6*CWBias);
% N=6;
% %  tau=0;
% ki=theta*v^2*N*tau.*(lambdaRdiff./(1+lambdaR./lambdaT)./(lambdaR*theta+Drot))./(1+tau.*(lambdaR+Drot))/3;
%  gradient=((f(3:end)-f(1:end-2))/2/dx);
%  plot((x(1:end-2)+x(3:end))/2,gradient);
%  hold on;
%  plot(x(phenotypemaxind),CWBias);
%  plot(x(phenotypemaxind),ki./max(ki));
%  phenotypemax_loc=smooth(smooth(x(phenotypemaxind)));
% %  phenotypemaxind=smooth(smooth(phenotypemaxind));
%  plot(phenotypemax_loc(1:7),ki(1:7),'.','MarkerEdgeColor',mycmap(1,:));
%  hold on;
%  plot(phenotypemax_loc(519:522),ki(519:522),'.','MarkerEdgeColor',mycmap(end,:));
%  for i=1:8
%     scatter(phenotypemax_loc(i+7:8:519),ki(i+7:8:519),25,mycmap);
%  end
%  %%
% %gradientO=((f(3:end)-f(1:end-2))/2/dx+(fO(3:end)-fO(1:end-2))/2/dx);
% plot(CWBias,gradient(phenotypemaxind)'.*ki);%+gradientO(phenotypemaxind)'.*ki/5000);
% % hold on;
% % profile1=sum(dyn_profile(i).cell_density(:,:,end));
% % profile2=sum(dyn_profile(i).cell_density(:,:,end-1));
% % [a, ind1]=max(profile1);
% % [a,ind2]=max(profile2);
% 
%  
% 
% plot(CWBias,ones(1,length(CWBias)).*(ind1-ind2)*dx/120,'--');
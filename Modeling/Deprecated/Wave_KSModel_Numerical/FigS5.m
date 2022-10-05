% load('C:\Users\XF\Dropbox\paper_wave_diversity\Wave_KSModel_Numerical\test4.mat');
% Aspind=log(20):(log(300)-log(20))/100:log(300);
% Timeind=log(20):(log(55)-log(20))/100:log(55);
% Asp=exp(Aspind);Time=15:40/(length(Asp)-1):55;tic;
% parfor i=1:length(Asp)
 Asp=[50 100 200];
 Time=[15 21 27];tic;
 for i=1:3
       [wavespeed_late(i) final_profile(i) dyn_profile(i) dx]=chemo_oxygen_fun_v2_no_oxg_grad(Asp(i),Time(i)*60);
 end
 %% peak position
for i=1:3
    cell_density=sum(dyn_profile(i).cell_density(:,:,:));
    clear celldensity
    celldensity(:,:)=cell_density(1,:,:);
    [a,b]=max(celldensity(50:end,:));
    t=0:0.2*50:0.2*(length(a)-1)*50;
    plot(smooth(b)*30/1000,t/60);hold on;
end
%%
for i=1:3
    cell_density=sum(dyn_profile(i).cell_density(:,:,:));
    clear celldensity
    celldensity(:,:)=cell_density(1,:,:);
    plot(celldensity(:,end)); hold on;
end



   
%%
% dx=20;
%  cell_density(i,:)=sum(final_profile(i).cell_density);
%     [denmax(i),denmaxind(i)]=max(cell_density(i,200:end));
%     wavedensity(i)=sum(cell_density(i,(200+denmaxind(i)-200):end));
%     cell_number(i)=wavedensity(i)*(dx*6*8.4/10*1.4);
%     
%   dyn_den_end=dyn_profile(i).cell_density(:,(200+denmaxind(i)-200):(200+denmaxind(i)+100),end);
%   [phenotypemax,phenotypemaxind]=max(dyn_den_end');
%   x=0:dx:(length((200+denmaxind(i)-200):(200+denmaxind(i)+100))-1)*dx;
%   plot(x,sum(dyn_den_end)*8.5*10^8./sum(sum(dyn_den_end)*8.5*10^8),'DisplayName','cell_density')
%   hold on;
%   dCWBias=0.001; CWBias=0.001:dCWBias:0.99;P = lognpdf(CWBias,log(0.295),0.308)/sum(lognpdf(CWBias,log(0.295),0.308))/dCWBias;
%   plot(x,dyn_den_end(abs(CWBias-0.1)<0.0009,:)*8.5*10^8/sum(dyn_den_end(abs(CWBias-0.1)<0.0009,:)*8.5*10^8));
%   plot(x,dyn_den_end(abs(CWBias-0.2)<0.0009,:)*8.5*10^8/sum(dyn_den_end(abs(CWBias-0.2)<0.0009,:)*8.5*10^8));
%   plot(x,dyn_den_end(abs(CWBias-0.3)<0.0009,:)*8.5*10^8/sum(dyn_den_end(abs(CWBias-0.3)<0.0009,:)*8.5*10^8));
%   plot(x,dyn_den_end(abs(CWBias-0.4)<0.0009,:)*8.5*10^8/sum(dyn_den_end(abs(CWBias-0.4)<0.0009,:)*8.5*10^8));
%  %%
% plot(x-x(100),sum(dyn_den_end)*8.5*10^8,'DisplayName','cell_density')
%   hold on;
%   dCWBias=0.001; CWBias=0.001:dCWBias:0.99;P = lognpdf(CWBias,log(0.295),0.308)/sum(lognpdf(CWBias,log(0.295),0.308))/dCWBias;
%   plot(x-x(100),sum(dyn_den_end(CWBias<0.2,:))*8.5*10^8);
%   plot(x-x(100),sum(dyn_den_end(abs(CWBias-0.25)<0.05,:))*8.5*10^8);
%   plot(x-x(100),sum(dyn_den_end(abs(CWBias-0.35)<0.05,:))*8.5*10^8);
%   plot(x-x(100),sum(dyn_den_end(abs(CWBias-0.45)<0.05,:))*8.5*10^8);
%  %%
% %   
% Ki=0.3;      %uM
% Ka=1000;    %uM
% KiO1=40;
% KiO2=150;
% KaO1=330;
% KaO2=10000;
%  x=0:dx:(length((200+denmaxind(i)-200):(200+denmaxind(i)+100))-1)*dx;
%  ASP_profile=dyn_profile(i).Asp((200+denmaxind(i)-200):(200+denmaxind(i)+100),end);
% %  O_profile=dyn_profile(i).O((200+denmaxind(i)-100):(200+denmaxind(i)+100),end);
%  f=log((1+ASP_profile/Ki)./(1+ASP_profile/Ka));
% %  fO=(log((1+O_profile/KiO1)./(1+O_profile/KaO1))-log((1+O_profile/KiO2)./(1+O_profile/KaO2)))/5000;
% %  plot((x(1:end-2)+x(3:end))/2,((f(3:end)-f(1:end-2))/2/dx+(fO(3:end)-fO(1:end-2))/2/dx)./max((f(3:end)-f(1:end-2))/2/dx+(fO(3:end)-fO(1:end-2))/2/dx)*5*10^9);
% plot(x-x(200),ASP_profile./200);hold on
%  plot((x(1:end-2)+x(3:end))/2-x(200),((f(3:end)-f(1:end-2))/2/dx)./max((f(3:end)-f(1:end-2))/2/dx));
%  %%
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
% 
%   dyn_den_end=dyn_profile(i).cell_density(:,(200+denmaxind(i)-100):(200+denmaxind(i)+100),end);
% 
%   x=0:dx:(length((200+denmaxind(i)-100):(200+denmaxind(i)+100))-1)*dx;
%   [phenotypemax,phenotypemaxind]=max(dyn_den_end');
%  ASP_profile=dyn_profile(i).Asp((200+denmaxind(i)-100):(200+denmaxind(i)+100),end);
%  f=log((1+ASP_profile/Ki)./(1+ASP_profile/Ka));
% ki=theta*v^2*N*tau.*(lambdaRdiff./(1+lambdaR./lambdaT)./(lambdaR*theta+Drot))./(1+tau.*(lambdaR+Drot))/3;
% gradient=smooth(((f(3:end)-f(1:end-2))/2/dx));%gradientO=((f(3:end)-f(1:end-2))/2/dx+(fO(3:end)-fO(1:end-2))/2/dx);
% plot(CWBias,smooth(gradient(phenotypemaxind)'.*ki));%+gradientO(phenotypemaxind)'.*ki/5000);
% hold on;
% profile1=sum(dyn_profile(i).cell_density(:,100:end,end));
% profile2=sum(dyn_profile(i).cell_density(:,100:end,end-1));
% [a, ind1]=max(profile1);
% [a,ind2]=max(profile2);
% plot(CWBias,ones(1,length(CWBias)).*(ind1-ind2)*dx/40,'--');
% % % % 
% % %%
% % %%%%%%%Figure 4A
% Ki=0.3;      %uM
% Ka=1000;    %uM
% den_profile=dyn_profile(end).cell_density(:,200:end,end-1);
% ASP_profile=dyn_profile(i).Asp(200:end,end-1);
% 
% [den_profile_max,den_profile_maxind]=max(sum(den_profile(:,:)));
% x=0:dx:(length(ASP_profile)-1)*dx;
% x=x-x(den_profile_maxind);
% figure;subplot(2,1,1);
% plot(x(den_profile_maxind-200:den_profile_maxind+100),ASP_profile(den_profile_maxind-200:den_profile_maxind+100));
%  f=log((1+ASP_profile/Ki)./(1+ASP_profile/Ka));
%  
% subplot(2,1,2);
% gradientf=(-f(den_profile_maxind-200:den_profile_maxind+100-2)+f(den_profile_maxind-200+2:den_profile_maxind+100))/2/dx;
% plot((x(den_profile_maxind-200:den_profile_maxind+100-2)+x(den_profile_maxind-200+2:den_profile_maxind+100))/2,smooth(gradientf));
% [gradientfmax,gradientfmaxind]=max(gradientf);
% hold on;
% plot(x(gradientfmaxind+den_profile_maxind-200),gradientfmax,'ro');
% plot(x(gradientfmaxind+den_profile_maxind-200)*ones(length(gradientf),1),gradientf,'r--');
%  %%
%  [phenotypemax,phenotypemaxind]=max(den_profile(:,den_profile_maxind-200:den_profile_maxind+100)');
%  X_max=x(phenotypemaxind+den_profile_maxind-200);
% plot(smooth(X_max(13:503)),ki(13:503));
% hold on;
% plot(ones(1,length(ki))*x(gradientfmaxind+den_profile_maxind-200),ki,'r--');
% hold on;
% plot(x,ki(503)*ones(1,length(x)),'r--');
% figure;
% plot(CWBias(13:503),X_max(13:503))
% %%
% gradientf=(-f(1:end-2)+f(3:end))/2/dx;
% gradient_F=gradientf(phenotypemaxind+den_profile_maxind-200-1);
% plot(smooth(X_max(13:503)),smooth(gradient_F(13:503).*ki(13:503)'));
% [a, ind1]=max(profile1);
% [a,ind2]=max(profile2);
% 
%  hold on;
% 
% plot(X_max(13:503),ones(1,length(X_max(13:503))).*(ind1-ind2)*dx/40,'--');
% ylim([0 6]);
% figure;
% plot(smooth(X_max(13:503)),smooth(gradient_F(13:503)));
% 
% %% x vs mean ki
% for i=1:length(x)
%     meanki(:,i)=sum(ki'.*den_profile(:,i)./sum(den_profile(:,:)')');
% end
% 
% %%
%  %%
% % %%%%%Figure 4B
% % [phenotypemax,phenotypemaxind]=max(den_profile(:,den_profile_maxind-200:den_profile_maxind+100)');
% % plot(CWBias,(x(phenotypemaxind+den_profile_maxind-200)));
% % hold on;
% % plot(CWBias,ones(1,length(CWBias))*x(gradientfmaxind+den_profile_maxind-100),'r--');
% % hold on;
% % plot(CWBias(gradientfmaxind+den_profile_maxind)*ones(1,length(x)),x,'r--');
% % for i=1:length(CWBias)
% %    phenotypemeanpos(i)=sum(den_profile(i,den_profile_maxind-200:end).*x(den_profile_maxind-200:end))./sum(den_profile(i,den_profile_maxind-200:end));
% % end
% % 
% % plot(phenotypemeanpos,CWBias);
% %  %%
% % %%%%%%Figure 4C
% % [gradmax,gradmaxind]=max(gradientf(phenotypemaxind));
% % figure;
% % subplot(2,1,1);
% % hold on;
% % plot(CWBias(1:gradmaxind),smooth(gradientf(phenotypemaxind(1:gradmaxind)))');
% % subplot(2,1,2);
% % plot(CWBias(1:gradmaxind),ki(1:gradmaxind));
% % figure;
% % plot(x(den_profile_maxind-100:den_profile_maxind+100)',sum(den_profile(:,den_profile_maxind-100:den_profile_maxind+100)));
% % % %% Figure 4D
% % figure;
% % 
% % profile1=sum(dyn_profile(end).cell_density(:,:,end));
% % profile2=sum(dyn_profile(end).cell_density(:,:,end-1));
% % [a, ind1]=max(profile1);
% % [a,ind2]=max(profile2);
% % plot(CWBias(1:gradmaxind),smooth(smooth(gradientf(phenotypemaxind((1:gradmaxind))))'.*ki((1:gradmaxind))));
% %  hold on;
% % 
% % plot(CWBias(1:gradmaxind),ones(1,length(CWBias(1:gradmaxind))).*(ind1-ind2)*dx/120,'--');
% % ylim([0 6]);
% % 

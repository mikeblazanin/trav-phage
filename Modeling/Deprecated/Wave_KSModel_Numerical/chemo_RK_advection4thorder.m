% clear all;
% dr=[19.8649697038349,59.5615356798170,99.2581016557992,138.954667631781,178.651233607763,218.347799583746,258.044365559728,297.740931535710,337.437497511692,377.134063487674,416.830629463656,456.527195439639,496.223761415621,535.920327391603,575.616893367585,615.313459343567,655.010025319549,694.706591295531,734.403157271514,774.099723247496,813.796289223478,853.492855199460,893.189421175442,932.885987151424,972.582553127407,1012.27911910339,1051.97568507937,1091.67225105535,1131.36881703134,1171.06538300732,1210.76194898330,1250.45851495928,1290.15508093526,1329.85164691125,1369.54821288723,1409.24477886321,1448.94134483919,1488.63791081517,1528.33447679116,1568.03104276714,1607.72760874312,1647.42417471910,1687.12074069509,1726.81730667107,1766.51387264705,1806.21043862303,1845.90700459901,1885.60357057500,1925.30013655098,1964.99670252696];
% P=[0.0668601905392834,0.0859603756959549,0.106609096872780,0.115683326454114,0.118739187872561,0.108081810524300,0.0888129011889620,0.0692441312063391,0.0518318347848702,0.0382581428177714,0.0269610627553021,0.0207251257991432,0.0155364691360261,0.0122145177839684,0.0102948273876606,0.00955730605453259,0.00799259632093305,0.00657189733185484,0.00583554050609490,0.00473061709833230,0.00370934413651132,0.00335940967239820,0.00267584384732532,0.00212231467835660,0.00171609569145213,0.00153792606413330,0.00124466429193423,0.00125689161929925,0.00102573690673202,0.000984979148848626,0.000894341658698417,0.000792253179428589,0.000977021681833297,0.000865811228180043,0.000848343617658590,0.000870081088529732,0.000810303043634091,0.000480165204778621,0.000450664351453499,0.000447753083033257,0.000366819820950523,0.000296367125180660,0.000323150794646889,0.000288021489042633,0.000300636985530349,0.000278123176413809,0.000245516970107096,0.000176811035389379,0.000119685479498848,3.39647982361595e-05];
% %load the original distribution
% for i=1:length(dr) %calculate the Yp value
%     options=optimset('TolX',10^-100);
%     Yp(i)=fzero(@(yp)findYp(yp,dr(i)),[0,6],options);
% end
% %Yp=[Yp(13),Yp(13)];dr=[dr(13),dr(13)];P=[0.5,0.5];
%% calculate the chemotactic coefficients based on the PLOS Comp. 
clear all;
%    dCWBias=0.1; CWBias=0.01:dCWBias:0.99;P = lognpdf(CWBias,log(0.28),0.3)/sum(lognpdf(CWBias,log(0.28),0.3))/dCWBias;
   CWBias=[0.18 0.9];P=[0.9 0.1];dCWBias=1;
dt=0.05;    %s
dx=20;      %um
% CWBias=[0.3 0.8]
%    
%  plot(CWBias,P)
%  P=[0.2 0.8];
g1=40;g2=40;
Kd      =  3.06;
w       =  1.3;
Yp=g2*Kd./(g2-g1/2+log(1./CWBias-1))-Kd;
deltaG  =  g1/4-g2/2*(Yp./(Kd+Yp));
lambdaT   =  w*exp(deltaG);
lambdaR  =  w*exp(-deltaG);
CWBias  =  ((lambdaR)./((lambdaT)+(lambdaR)));
alpha=6;
a=Yp/alpha;
Drot=0.062*2;
tauR=1./(Drot+lambdaR*(1-0.0));
v=30;
Diff=v^2*(1-CWBias).*tauR/3;
lambdaRdiff=-(lambdaR).*a.*(a-1).*(alpha*g2*Kd)./(2*(Kd+alpha*a).^2); 
tau=20.8*exp(-5.86*CWBias);
N=6;
%  tau=0;
ki=v^2*N*tau.*(lambdaRdiff./(1+lambdaR./lambdaT)./(lambdaR*(1-0.0)+Drot))./(1+tau.*(lambdaR+Drot))/3*(1-0.0);
%%   

Ki=1.8;      %uM
Ka=500;    %uM
Asp=50;   %initial Asp conc in units of uM 
L=14000/dx;
LX=0:dx:(L*dx-dx);
u=zeros(round(L),length(CWBias));
Drho=ones(round(L),1)*Diff/dx^2;
A=ones(round(L),1)*Asp;

T=3100;%
Dn=500/dx^2;
G0=1/6;% 10uM/min for OD 1 cells.
DL=1500/dx; V=250/dx;
for x=1:round(DL)
    ux(x)=exp(-(x)^2/V^2);
%     Ax(x)=(x/(V*2)).^5./(1+(x/(V*2)).^5);
end

u(3:round(DL)+2,:)=ux'*P*75*dCWBias; %unit is OD1,~5*dx cells in 600um*10um*dx %total cell: 0.7*10^9*600*10*10^4*2*10^-12=8*10^4 cells
u(1,:)=u(3,:);u(2,:)=u(3,:);
ua=u;ub=u;uc=u;ud=u;Aa=A;Ab=A;Ac=A;Ad=A;
% 
% A(3:round(DL)+2)=Asp*Ax;
% A(1)=A(3);A(2)=A(3);A(end)=A(end-3);A(end-1)=A(end-2);
% figure;
% plot(LX,sum(u')',LX,A);
sum(sum(u))
onemaxtrix=ones(1,length(CWBias));
%%

tic
% for t=1:100/dt
%      u1=u;A1=A;f=u;
%       Gamma1=G0*0*sum(u1')'.*(0+1./(1+(sum(u1')'/10).^3))./(1+1./A1);
% %    f=log((1+A1/Ki)./(1+A1/Ka))*ki/dx^2;
% %     ua(2:end-1,:)=-u1(2:end-1,:).*(2*Drho(2:end-1,:))+Drho(2:end-1,:).*(u1(1:end-2,:)+u1(3:end,:))-(u1(3:end,:)-u1(1:end-2,:)).*((f(3:end,:)-f(1:end-2,:)))/4-f(2:end-1,:).*((A1(3:end)+A1(1:end-2)-2*A1(2:end-1))*onemaxtrix);
%     Aa(2:end-1)=-A1(2:end-1)*(2*Dn)+Dn*(A1(1:end-2)+A1(3:end))-Gamma1(2:end-1);
%    Aa(1)=Aa(4);Aa(2)=Aa(3);Aa(end)=Aa(end-3);Aa(end-1)=Aa(end-2);
%     A1=A+Aa*dt/2;
%     
% %%
%     Gamma1=G0*0*sum(u1')'.*(0+1./(1+(sum(u1')'/10).^3))./(1+1./A1);
%      Ab(2:end-1)=-A1(2:end-1)*(2*Dn)+Dn*(A1(1:end-2)+A1(3:end))-Gamma1(2:end-1);
%    Ab(1)=Ab(4);Ab(2)=Ab(3);Ab(end)=Ab(end-3);Ab(end-1)=Ab(end-2);
%     Ab(end)=Ab(end-1);
%     %%
%     
%      Gamma1=G0*0*sum(u1')'.*(0+1./(1+(sum(u1')'/10).^3))./(1+1./A1);
%       Ac(2:end-1)=-A1(2:end-1)*(2*Dn)+Dn*(A1(1:end-2)+A1(3:end))-Gamma1(2:end-1);
%      Ac(1)=Ac(4);Ac(2)=Ac(3);Ac(end)=Ac(end-3);Ac(end-1)=Ac(end-2);
%     
%     %%
%   A1=A+Ac*dt;
%      Gamma1=G0*0*sum(u1')'.*(0+1./(1+(sum(u1')'/10).^3))./(1+1./A1);
%     Ad(2:end-1)=-A1(2:end-1)*(2*Dn)+Dn*(A1(1:end-2)+A1(3:end))-Gamma1(2:end-1);
%     Ad(1)=Ad(4);Ad(2)=Ad(3);Ad(end)=Ad(end-3);Ad(end-1)=Ad(end-2);
% %     %
% %     
%      
%      A(2:end-1)=A(2:end-1)+1/6*(Aa(2:end-1)+2*Ab(2:end-1)+2*Ac(2:end-1)+Ad(2:end-1))*dt;
%     A(1)=A(4);A(2)=A(3);A(end)=A(end-3);A(end-1)=A(end-2);
%     if (mod(t,60/dt)==0)
% %          plot(LX,u1(:,1),'g');hold on;plot(Lx,u1(:,2),'b');
% %         u(:,2)=0;
% %         plot(LX,sum(u')','g');hold on;%plot(n,'g');
%        plot(A,'r');hold on;
%         end
% end
%%
 figure;

for t=1:T
    u1=u;A1=A;
    Gamma1=G0*sum(u1')'.*(0+1./(1+(sum(u1')'/10).^3))./(1+1./A1);
    %Gamma1=G0*sum(u1')'./(1+0.5./A1);
    f=log((1+A1/Ki)./(1+A1/Ka))*ki/dx^2;
    g()=
    ua(3:end-2,:)=-u1(3:end-2,:).*(2*Drho(3:end-2,:))+Drho(3:end-2,:).*(u1(2:end-3,:)+u1(4:end-1,:))-(u1(4:end-1,:)-u1(2:end-3,:)).*(-f(5:end,:)+8*f(4:end-1,:)-8*f(2:end-3,:)+f(1:end-4,:))/12/2*0-u1(3:end-2,:).*(f(4:end-1,:)+f(2:end-3,:)-2*f(3:end-2,:));
    Aa(2:end-1)=-A1(2:end-1)*(2*Dn)+Dn*(A1(1:end-2)+A1(3:end))-Gamma1(2:end-1);
    Aa(1)=Aa(5);Aa(2)=Aa(4);Aa(3)=(Aa(2)+Aa(4))/2;
    Aa(end)=Aa(end-3);Aa(end-1)=Aa(end-2);
    ua(1,:)=ua(5,:);ua(end,:)=ua(end-3,:);ua(2,:)=ua(4,:);ua(end-1,:)=ua(end-2,:);ua(3,:)=(ua(2,:)+ua(4,:))/2;
    %%
    u1=u+ua*dt/2;A1=A+Aa*dt/2;
    Gamma1=G0*sum(u1')'.*(0+1./(1+(sum(u1')'/10).^3))./(1+1./A1);
     f=log((1+A1/Ki)./(1+A1/Ka))*ki/dx^2;
    ub(3:end-2,:)=-u1(3:end-2,:).*(2*Drho(3:end-2,:))+Drho(3:end-2,:).*(u1(2:end-3,:)+u1(4:end-1,:))-(u1(4:end-1,:)-u1(2:end-3,:)).*(-f(5:end,:)+8*f(4:end-1,:)-8*f(2:end-3,:)+f(1:end-4,:))/12/2*0-u1(3:end-2,:).*(f(4:end-1,:)+f(2:end-3,:)-2*f(3:end-2,:));
    Ab(2:end-1)=-A1(2:end-1)*(2*Dn)+Dn*(A1(1:end-2)+A1(3:end))-Gamma1(2:end-1);
    Ab(1)=Ab(5);Ab(2)=Ab(4);Ab(3)=(Ab(2)+Ab(4))/2;
    Ab(end)=Ab(end-3);Ab(end-1)=Ab(end-2);
    ub(1,:)=ub(5,:);    ub(end,:)=ub(end-3,:);ub(2,:)=ub(4,:);ub(end-1,:)=ub(end-2,:);ub(3,:)=(ub(2,:)+ub(4,:))/2;
    %%
    u1=u+ub*dt/2;A1=A+Ab*dt/2;
    Gamma1=G0*sum(u1')'.*(0+1./(1+(sum(u1')'/10).^3))./(1+1./A1);
     f=log((1+A1/Ki)./(1+A1/Ka))*ki/dx^2;
    uc(3:end-2,:)=-u1(3:end-2,:).*(2*Drho(3:end-2,:))+Drho(3:end-2,:).*(u1(2:end-3,:)+u1(4:end-1,:))-(u1(4:end-1,:)-u1(2:end-3,:)).*(-f(5:end,:)+8*f(4:end-1,:)-8*f(2:end-3,:)+f(1:end-4,:))/12/2*0-u1(3:end-2,:).*(f(4:end-1,:)+f(2:end-3,:)-2*f(3:end-2,:));
    Ac(2:end-1)=-A1(2:end-1)*(2*Dn)+Dn*(A1(1:end-2)+A1(3:end))-Gamma1(2:end-1);
     Ac(1)=Ac(5);Ac(2)=Ac(4);Ac(3)=(Ac(2)+Ac(4))/2;
     Ac(end)=Ac(end-3);Ac(end-1)=Ac(end-2);
    uc(1,:)=uc(5,:);uc(end,:)=uc(end-3,:);uc(2,:)=uc(4,:);uc(end-1,:)=uc(end-2,:);uc(3,:)=(uc(2,:)+uc(4,:))/2;
    
    %%
    u1=u+uc*dt;A1=A+Ac*dt;
    Gamma1=G0*sum(u1')'.*(0+1./(1+(sum(u1')'/10).^3))./(1+1./A1);
     f=log((1+A1/Ki)./(1+A1/Ka))*ki/dx^2;
    ud(3:end-2,:)=-u1(3:end-2,:).*(2*Drho(3:end-2,:))+Drho(3:end-2,:).*(u1(2:end-3,:)+u1(4:end-1,:))-(u1(4:end-1,:)-u1(2:end-3,:)).*(-f(5:end,:)+8*f(4:end-1,:)-8*f(2:end-3,:)+f(1:end-4,:))/12/2*0-u1(3:end-2,:).*(f(4:end-1,:)+f(2:end-3,:)-2*f(3:end-2,:));
    Ad(2:end-1)=-A1(2:end-1)*(2*Dn)+Dn*(A1(1:end-2)+A1(3:end))-Gamma1(2:end-1);
     Ad(1)=Ad(5);Ad(2)=Ad(4);Ad(3)=(Ad(2)+Ad(4))/2;
     Ad(end)=Ad(end-3);Ad(end-1)=Ad(end-2);
    ud(1,:)=ud(5,:);ud(end,:)=ud(end-3,:);ud(2,:)=ud(4,:);ud(end-1,:)=ud(end-2,:);ud(3,:)=(ud(2,:)+ud(4,:))/2;
%     %
%     
     u=u+1/6*(ua+2*ub+2*uc+ud)*dt;
     A=A+1/6*(Aa+2*Ab+2*Ac+Ad)*dt;
     A(1)=A(5);A(2)=A(4);A(3)=(A(2)+A(4))/2;
     A(end)=A(end-3);A(end-1)=A(end-2);

     u(1,:)=u(5,:);u(end,:)=u(end-3,:);u(2,:)=u(4,:);u(end-1,:)=u(end-2,:);u(3,:)=(u(2,:)+u(4,:))/2;
     [maxden(t),ind(t)]=max(sum(u(4:end-2,:)'));
     tot(t)=sum(sum(u(4:end-2,:)));    
    %%    
%     u1=u;A1=A;
%      Gamma1=G0*sum(u1')'.*(0+1./(1+(sum(u1')'/10).^3))./(1+1./A1);
%     %Gamma1=G0*sum(u1')'./(1+0.5./A1);
%      f=log((1+A1/Ki)./(1+A1/Ka))*ki/dx^2;
%     ua(2:end-1,:)=-u1(2:end-1,:).*(2*Drho(2:end-1,:))+Drho(2:end-1,:).*(u1(1:end-2,:)+u1(3:end,:))-(u1(3:end,:)-u1(1:end-2,:)).*((f(3:end,:)-f(1:end-2,:)))/4-f(2:end-1,:).*(f(3:end,:)+f(1:end-2,:)-2*f(2:end-1,:));
%     Aa(2:end-1)=-A1(2:end-1)*(2*Dn)+Dn*(A1(1:end-2)+A1(3:end))-Gamma1(2:end-1);
%     Aa(1)=Aa(2);ua(1,:)=ua(2,:);Aa(end)=Aa(end-1);ua(end,:)=ua(end-1,:);
%     %%
%     u1=u+ua*dt/2;A1=A+Aa*dt/2;
%      Gamma1=G0*sum(u1')'.*(0+1./(1+(sum(u1')'/10).^3))./(1+1./A1);
%       f=log((1+A1/Ki)./(1+A1/Ka))*ki/dx^2;
%     ub(2:end-1,:)=-u1(2:end-1,:).*(2*Drho(2:end-1,:))+Drho(2:end-1,:).*(u1(1:end-2,:)+u1(3:end,:))-(u1(3:end,:)-u1(1:end-2,:)).*((f(3:end,:)-f(1:end-2,:)))/4-f(2:end-1,:).*(f(3:end,:)+f(1:end-2,:)-2*f(2:end-1,:));
%     Ab(2:end-1)=-A1(2:end-1)*(2*Dn)+Dn*(A1(1:end-2)+A1(3:end))-Gamma1(2:end-1);
%     Ab(1)=Ab(2);ub(1,:)=ub(2,:);Ab(end)=Ab(end-1);ub(end,:)=ub(end-1,:);
%     %%
%     u1=u+ub*dt/2;A1=A+Ab*dt/2;
%       Gamma1=G0*sum(u1')'.*(0+1./(1+(sum(u1')'/10).^3))./(1+1./A1);
%     f=log((1+A1/Ki)./(1+A1/Ka))*ki/dx^2;
%     uc(2:end-1,:)=-u1(2:end-1,:).*(2*Drho(2:end-1,:))+Drho(2:end-1,:).*(u1(1:end-2,:)+u1(3:end,:))-(u1(3:end,:)-u1(1:end-2,:)).*((f(3:end,:)-f(1:end-2,:)))/4-f(2:end-1,:).*(f(3:end,:)+f(1:end-2,:)-2*f(2:end-1,:));
%     Ac(2:end-1)=-A1(2:end-1)*(2*Dn)+Dn*(A1(1:end-2)+A1(3:end))-Gamma1(2:end-1);
%     Ac(1)=Ac(2);uc(1,:)=uc(2,:);Ac(end)=Ac(end-1);uc(end,:)=uc(end-1,:);
%     
%     %%
%     u1=u+uc*dt;A1=A+Ac*dt;
%       Gamma1=G0*sum(u1')'.*(0+1./(1+(sum(u1')'/10).^3))./(1+1./A1);
%      f=log((1+A1/Ki)./(1+A1/Ka))*ki/dx^2;
%     ud(2:end-1,:)=-u1(2:end-1,:).*(2*Drho(2:end-1,:))+Drho(2:end-1,:).*(u1(1:end-2,:)+u1(3:end,:))-(u1(3:end,:)-u1(1:end-2,:)).*((f(3:end,:)-f(1:end-2,:)))/4-f(2:end-1,:).*(f(3:end,:)+f(1:end-2,:)-2*f(2:end-1,:));
%     Ad(2:end-1)=-A1(2:end-1)*(2*Dn)+Dn*(A1(1:end-2)+A1(3:end))-Gamma1(2:end-1);
%     Ad(1)=Ad(2);ud(1,:)=ud(2,:);
%     Ad(end)=Ad(end-1);ud(end,:)=ud(end-1,:);
% %     %
% %     
%      u(2:end-1,:)=u(2:end-1,:)+1/6*(ua(2:end-1,:)+2*ub(2:end-1,:)+2*uc(2:end-1,:)+ud(2:end-1,:))*dt;
%      A(2:end-1)=A(2:end-1)+1/6*(Aa(2:end-1)+2*Ab(2:end-1)+2*Ac(2:end-1)+Ad(2:end-1))*dt;
%           A(1)=A(2);u(1,:)=u(2,:);
%      A(end)=A(end-1);u(end,:)=u(end-1,:);
     %%
%     u(2:end-1,:)=u1(2:end-1,:).*(1-2*Drho(2:end-1,:)*dt)+dt*Drho(2:end-1,:).*(u1(1:end-2,:)+u1(3:end,:))-dt*(u1(3:end,:)-u1(1:end-2,:)).*((f(3:end,:)-f(1:end-2,:)))/4-dt*f(2:end-1,:).*((A1(3:end)+A1(1:end-2)-2*A1(2:end-1))*onemaxtrix);
%     A(2:end-1)=A1(2:end-1)*(1-2*Dn*dt)+dt*Dn*(A1(1:end-2)+A1(3:end))-dt*Gamma1(2:end-1);
%     f(2:end-1,:)=0* 0*(Ka-Ki)*u1(2:end-1,:).*((1./((Ki+A1(2:end-1)).*(Ka+A1(2:end-1))))*ki)/dx.*((f(3:end,:)-f(1:end-2,:)))/2;
%     f(2:end-1,:)=u1(2:end-1,:);
%     f(1,:)=u1(1,:);f(end,:)=u1(end,:);
%     utemp=(u(1:end-1,:)+u(2:end,:))/2-dt/2*(f(2:end,:)-f(1:end-1,:))/dx;
%     ftemp=utemp;  
%     fplus(2:end-1,:)=0* 0*(Ka-Ki)*u1.*((1./((Ki+A1(2:end-1)).*(Ka+A1(2:end-1))))*ki)/dx.*((f(3:end,:)-f(1:end-2,:)))/2;fplus(1,:)=0;fplus(end,:)=0;
  
%     fminus(2:end-1,:)=0* 0*(Ka-Ki)*u1.*((1./((Ki+A1(2:end-1)).*(Ka+A1(2:end-1))))*ki)/dx.*((f(3:end,:)-f(1:end-2,:)))/2;fminus(1,:)=0;fminus(end,:)=0;
%     u(2:end-1,:)=u(2:end-1,:).*(1-2*Drho(2:end-1,:)*dt)+dt*Drho(2:end-1,:).*(u(1:end-2,:)+u(3:end,:))-dt*(ftemp(2:end,:)-ftemp(1:end-1,:))/dx;%dt*(u1(3:end,:)-u1(1:end-2,:)).*((f(3:end,:)-f(1:end-2,:)))/4-dt*f(2:end-1,:).*((A1(3:end)+A1(1:end-2)-2*A1(2:end-1))*onemaxtrix);
%     A(2:end-1)=A1(2:end-1)*(1-2*Dn*dt)+dt*Dn*(A1(1:end-2)+A1(3:end))-dt*Gamma1(2:end-1);

%     A(1)=A(2);u(1,:)=u(2,:);
%     A(end)=A(end-1);u(end,:)=u(end-1,:);
%       tot(t)=sum(sum(u(2:end-1,:))); 
%     if t==3400
%         plot(LX,sum(u')'/0.84,'g');hold on; plot(LX,A,'r');
%     end
%     if (mod(t,600/dt)==0)
%          plot(LX,u1(:,1),'g');hold on;plot(Lx,u1(:,2),'b');
%         u(:,2)=0;
%         plot(LX,sum(u')'/0.84,'g');hold on;%plot(n,'g');
%       plot(A,'y');
%         toc
%      end
end
% figure;plot(CWBias,sum(u(500:1000,:))/sum(sum(u(500:1000,:)))/0.01)
% hold on;
% plot(CWBias, P)
    %c=(sum(ur>0.4)).2;
%c=c';
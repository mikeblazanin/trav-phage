% clear all;
% dr=[19.8649697038349,59.5615356798170,99.2581016557992,138.954667631781,178.651233607763,218.347799583746,258.044365559728,297.740931535710,337.437497511692,377.134063487674,416.830629463656,456.527195439639,496.223761415621,535.920327391603,575.616893367585,615.313459343567,655.010025319549,694.706591295531,734.403157271514,774.099723247496,813.796289223478,853.492855199460,893.189421175442,932.885987151424,972.582553127407,1012.27911910339,1051.97568507937,1091.67225105535,1131.36881703134,1171.06538300732,1210.76194898330,1250.45851495928,1290.15508093526,1329.85164691125,1369.54821288723,1409.24477886321,1448.94134483919,1488.63791081517,1528.33447679116,1568.03104276714,1607.72760874312,1647.42417471910,1687.12074069509,1726.81730667107,1766.51387264705,1806.21043862303,1845.90700459901,1885.60357057500,1925.30013655098,1964.99670252696];
% P=[0.0668601905392834,0.0859603756959549,0.106609096872780,0.115683326454114,0.118739187872561,0.108081810524300,0.0888129011889620,0.0692441312063391,0.0518318347848702,0.0382581428177714,0.0269610627553021,0.0207251257991432,0.00155364691360261,0.00122145177839684,0.00102948273876606,0.00955730605453259,0.00799259632093305,0.00657189733185484,0.00583554050609490,0.00473061709833230,0.00370934413651132,0.00335940967239820,0.00267584384732532,0.00212231467835660,0.00171609569145213,0.00153792606413330,0.00124466429193423,0.00125689161929925,0.00102573690673202,0.001984979148848626,0.001894341658698417,0.001792253179428589,0.001977021681833297,0.001865811228180043,0.001848343617658590,0.001870081088529732,0.001810303043634091,0.001480165204778621,0.001450664351453499,0.001447753083033257,0.001366819820950523,0.001296367125180660,0.001323150794646889,0.001288021489042633,0.001300636985530349,0.001278123176413809,0.001245516970107096,0.00176811035389379,0.00119685479498848,3.39647982361595e-05];
% %load the original distribution
% for i=1:length(dr) %calculate the Yp value
%     options=optimset('TolX',10^-100);
%     Yp(i)=fzero(@(yp)findYp(yp,dr(i)),[0,6],options);
% end
% %Yp=[Yp(13),Yp(13)];dr=[dr(13),dr(13)];P=[0.5,0.5];
%% calculate the chemotactic coefficients based on the PLOS Comp. 
clear all;
% dCWBias=0.1; CWBias=0.001:dCWBias:0.99;P = lognpdf(CWBias,log(0.28),0.3)/sum(lognpdf(CWBias,log(0.28),0.3))/dCWBias;
    CWBias=[0.1;0.1];P=[0.9;0.1];dCWBias=1;
    dt=0.02;    %s
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
tauR=1./(Drot+lambdaR);
v=30;
Diff=v^2*(1-CWBias).*tauR/3;
lambdaRdiff=-(lambdaR).*a.*(a-1).*(alpha*g2*Kd)./(2*(Kd+alpha*a).^2); 
tau=20.8*exp(-5.86*CWBias);
N=6;
%  tau=0;
ki=v^2*N*tau.*(lambdaRdiff./(1+lambdaR./lambdaT)./(lambdaR*(1-0.1)+Drot))./(1+tau.*(lambdaR+Drot))/3;
%%   

Ki=1.8;      %uM
Ka=1000;    %uM
Asp=200;   %initial Asp conc in units of uM 
L=40000/dx;
LX=0:dx:(L*dx-dx);
u=zeros(length(CWBias),round(L));
Drho=Diff/dx^2*ones(1,round(L));
A=ones(1,round(L))*Asp;

T=60*40/dt;%
Dn=800/dx^2;
G0=1/6/0.8;% 10uM/min for OD 1 cells.
DL=1500/dx; V=250/dx;
for x=1:round(DL)
    ux(x)=exp(-(x)^2/V^2);
%      Ax(x)=(x/(V*2)).^5./(1+(x/(V*2)).^5);
end

u(:,2:round(DL)+1)=P*ux*80*dCWBias; %unit is OD1,~5*dx cells in 600um*10um*dx %total cell: 0.7*10^9*600*10*10^4*2*10^-12=8*10^4 cells
u(:,1)=u(:,2);
ua=u;ub=u;uc=u;ud=u;Aa=A;Ab=A;Ac=A;Ad=A;
% 
%  A(2:round(DL)+1)=Asp*Ax;A(1)=A(2);
% A(1)=A(3);A(2)=A(3);A(end)=A(end-3);A(end-1)=A(end-2);
% figure;
plot(LX,u,LX,A);
sum(sum(u(:,2:end-1)))

tic
%%
g=u;grad=A;
for t=1:T
    u1=u;
    %4th RK for A;
        A1=A;
        Gamma1=G0*sum(u1);
        Aa(2:end-1)=Dn*(A1(1:end-2)+A1(3:end)-2*A1(2:end-1))-Gamma1(2:end-1);
        Aa(1)=Aa(2); Aa(end)=Aa(end-1); 
        A1=A+Aa*dt/2;A1(A1<=0)=0;
        Gamma1=G0*sum(u1);
        Ab(2:end-1)=Dn*(A1(1:end-2)+A1(3:end)-2*A1(2:end-1))-Gamma1(2:end-1);
        Ab(1)=Ab(2); Ab(end)=Ab(end-1); 
        A1=A+Ab*dt/2;A1(A1<=0)=0;
        Gamma1=G0*sum(u1);
        Ac(2:end-1)=Dn*(A1(1:end-2)+A1(3:end)-2*A1(2:end-1))-Gamma1(2:end-1);
        Ac(1)=Ac(2);Ac(end)=Ac(end-1);  
        A1=A+Ac*dt;A1(A1<=0)=0;
        Gamma1=G0*sum(u1);
        Ad(2:end-1)=Dn*(A1(1:end-2)+A1(3:end)-2*A1(2:end-1))-Gamma1(2:end-1);
        Ad(1)=Ad(2); Ad(end)=Ad(end-1); 
        A=A+1/6*(Aa+2*Ab+2*Ac+Ad)*dt;
        A(1)=A(2); A(end)=A(end-1);  
        A(A<=0)=0;
     % 4th RK for u;
     u1=u;A1=A;grad(2:end-1)=(A1(3:end)-A1(1:end-2))./A1(2:end-1)/2/dx; grad(isnan(grad))=0;grad(A1<=0.05)=0;grad(1)=grad(2);grad(end)=grad(end-1);
     g(:,2:end-1)=u1(:,2:end-1).*(ki*grad(2:end-1)); g(:,1)=-g(:,2);g(:,end)=-g(:,end-1);
     ua(:,2:end-1)=Drho(:,2:end-1).*(u1(:,1:end-2)+u1(:,3:end)-2*u1(:,2:end-1))-(g(:,3:end)-g(:,1:end-2))/2/dx; ua(:,1)=ua(:,2);ua(:,end)=ua(:,end-1);
     u1=u+ua*dt/2;
     g(:,2:end-1)=u1(:,2:end-1).*(ki*grad(2:end-1)); g(:,1)=-g(:,2);g(:,end)=-g(:,end-1);
     ub(:,2:end-1)=Drho(:,2:end-1).*(u1(:,1:end-2)+u1(:,3:end)-2*u1(:,2:end-1))-(g(:,3:end)-g(:,1:end-2))/2/dx; ub(:,1)=ub(:,2);ub(:,end)=ub(:,end-1);
     u1=u+ub*dt/2;
     g(:,2:end-1)=u1(:,2:end-1).*(ki*grad(2:end-1)); g(:,1)=-g(:,2);g(:,end)=-g(:,end-1);
     uc(:,2:end-1)=Drho(:,2:end-1).*(u1(:,1:end-2)+u1(:,3:end)-2*u1(:,2:end-1))-(g(:,3:end)-g(:,1:end-2))/2/dx; uc(:,1)=uc(:,2);uc(:,end)=uc(:,end-1);
     u1=u+uc*dt;
     g(:,2:end-1)=u1(:,2:end-1).*(ki*grad(2:end-1)); g(:,1)=-g(:,2);g(:,end)=-g(:,end-1);
     ud(:,2:end-1)=Drho(:,2:end-1).*(u1(:,1:end-2)+u1(:,3:end)-2*u1(:,2:end-1))-(g(:,3:end)-g(:,1:end-2))/2/dx;ud(:,1)=ud(:,2);ud(:,end)=ud(:,end-1);
     u=u+1/6*(ua+2*ub+2*uc+ud)*dt;
     u(:,1)=u(:,2);u(:,end)=u(:,end-1);
     sum(u);
     indx=find(diff(sign(diff(sum(u'))))==-2)+1;
     if ~isempty(indx)
        peakpos(t)=LX(indx(end));
     else 
         peakpos=LX(1);
     end
      tot(t)=sum(sum(u(:,2:end-1)));   
    if (mod(t,300/dt)==0)
%          plot(LX,u1(:,1),'g');hold on;plot(LX,u1(:,2),'b');
      
        plot(LX,u/0.84,'g');hold on;%plot(n,'g');
        plot(LX,A/20,'r'); 
        toc
     end
 
end
sum(sum(u(:,2:end-1)))

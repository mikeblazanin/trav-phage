%% calculate the chemotactic coefficients based on the PLOS Comp. 
 function [wavespeed final_profile dyn_profile dx]=chemo_oxygen_fun_v2_finegrid_no_oxg_grad(Asp,Time)
 dCWBias=0.005; CWBias=0.01:dCWBias:0.99;P = lognpdf(CWBias,log(0.295),0.308)/sum(lognpdf(CWBias,log(0.295),0.308))/dCWBias;

%      CWBias=[0.2 0.2];P=[0.5 0.5];dCWBias=1;
CC=winter(length(CWBias));
    dt=0.01;    %s
dx=8;      %um
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
v=25*4/pi;%2D run speed 26um/s
theta=1-0.33;
Diff=v^2*(1-CWBias).*tauR/3/theta;
lambdaRdiff=-(lambdaR).*a.*(a-1).*(alpha*g2*Kd)./(2*(Kd+alpha*a).^2); 
tau=20*exp(-6*CWBias);
N=6;
%  tau=0;
ki=theta*v^2*N*tau.*(lambdaRdiff./(1+lambdaR./lambdaT)./(lambdaR*theta+Drot))./(1+tau.*(lambdaR+Drot))/3;
kiO=ki/5000;
%%   

Ki=0.3;      %uM
Ka=1000;    %uM
KiO1=40;
KiO2=150;
KaO1=330;
KaO2=10000;% Studies of bacterial aerotaxis in a microfluidic device 15% transition from attractant to reppellent
%  Asp=100;   %initial Asp conc in units of uM
kaO=0.8;
ktO=80;
L=13000/dx;
LX=(0:dx:(L*dx-dx))/1000;
u=zeros(length(CWBias),round(L));
Drho=Diff'/dx^2*ones(1,round(L));
A=ones(1,round(L))*Asp;
Oex=250;
O=ones(1,round(L))*Oex;
KO=0.02;
GOx=8*4/60;
  T=round(Time/dt);%
Dn=600/dx^2;%A method for the determination of diffusion coefficients for small molecules in aqueous solution, Analytical Biochemistry,Volume 166, Issue 2, 1 November 1987, Pages 335-341
DO=2500/dx^2;%http://compost.css.cornell.edu/oxygen/oxygen.diff.water.html
G0=9.0*0.84/60;% 10uM/min for OD 1 cells.
DL=900/dx; V=900/dx;
for x=1:round(DL)
    ux(x)=1;%exp(-(x)^2/V^2);
end
ux=ux./sum(ux);
u(:,2:round(DL)+1)=P'*ux*500*50/dx*dCWBias; %unit is OD1,~6*8.4*1.4*dx cells in 600um*14um*dx %total cell: 0.7*10^9*(600*14*10^4+900*14*900+?)*10^-12=15~20*10^4 cells
u(:,1)=u(:,2);
A=Asp*(1-1./(1+(LX*1000/(V*dx*3)).^3));
O=Oex*(1-0.8./(1+(LX*1000/(V*dx*2)).^6));
ua=u;ub=u;uc=u;ud=u;Aa=A;Ab=A;Ac=A;Ad=A;Oa=O;Ob=O;Oc=O;Od=O;


%%
tt=1;
g=u;grad=A;gO=u;
for t=1:T
    u1=u;A1=A;A1(A1<=0)=0;
    O1=O; O1(O1<=0)=0;
%     Gamma1=G0*u1.*(0.1+0.9./(1+(u1/9.5).^3)).*(A1>0.001);
    Gamma1=G0*sum(u1)./(1+(0.5./A1)).*((kaO./(1+((Oex-O1)/ktO).^5))+1-kaO);
    GammaO=GOx*sum(u1)./(1+10./O1);


    %f=ki'*A1;
     f=ki'*log((1+A1/Ki)./(1+A1/Ka));%
    % fO=kiO'*(log((1+O1/KiO1)./(1+O1/KaO1))-log((1+O1/KiO2)./(1+O1/KaO2)));
    g(:,2:end-1)=u1(:,2:end-1).*(f(:,3:end)-f(:,1:end-2))/2/dx;   g(:,1)=-g(:,2);g(:,end)=-g(:,end-1);
%      gO(:,2:end-1)=u1(:,2:end-1).*(fO(:,3:end)-fO(:,1:end-2))/2/dx;   gO(:,1)=-gO(:,2);gO(:,end)=-gO(:,end-1);
    ua(:,2:end-1)=Drho(:,2:end-1).*(u1(:,1:end-2)+u1(:,3:end)-2*u1(:,2:end-1))-(g(:,3:end)-g(:,1:end-2))/2/dx;%-(gO(:,3:end)-gO(:,1:end-2))/2/dx;
    Aa(2:end-1)=Dn*(A1(1:end-2)+A1(3:end)-2*A1(2:end-1))-Gamma1(2:end-1);
    Oa(2:end-1)=DO*(O1(1:end-2)+O1(3:end)-2*O1(2:end-1))-GammaO(2:end-1)+KO*((Oex-O1(2:end-1))+abs(Oex-O1(2:end-1)))/2;
    Aa(1)=Aa(2); Aa(end)=Aa(end-1);    ua(:,1)=ua(:,2);ua(:,end)=ua(:,end-1);  Oa(1)=Oa(2); Oa(end)=Oa(end-1);   
 
    %%
    u1=u+ua*dt/2;A1=A+Aa*dt/2;A1(A1<=0)=0;
    Gamma1=G0*sum(u1)./(1+(0.5./A1)).*((kaO./(1+((Oex-O1)/ktO).^5))+1-kaO);
    O1=O+Oa*dt/2;O1(O1<=0)=0;
    GammaO=GOx*sum(u1)./(1+10./O1);


    %f=ki'*A1;
    f=ki'*log((1+A1/Ki)./(1+A1/Ka));%
    % fO=kiO'*(log((1+O1/KiO1)./(1+O1/KaO1))-log((1+O1/KiO2)./(1+O1/KaO2)));
    g(:,2:end-1)=u1(:,2:end-1).*(f(:,3:end)-f(:,1:end-2))/2/dx;   g(:,1)=-g(:,2);g(:,end)=-g(:,end-1);
%      gO(:,2:end-1)=u1(:,2:end-1).*(fO(:,3:end)-fO(:,1:end-2))/2/dx;   gO(:,1)=-gO(:,2);gO(:,end)=-gO(:,end-1);
    ub(:,2:end-1)=Drho(:,2:end-1).*(u1(:,1:end-2)+u1(:,3:end)-2*u1(:,2:end-1))-(g(:,3:end)-g(:,1:end-2))/2/dx;%-(gO(:,3:end)-gO(:,1:end-2))/2/dx;
    Ab(2:end-1)=Dn*(A1(1:end-2)+A1(3:end)-2*A1(2:end-1))-Gamma1(2:end-1);
    Ob(2:end-1)=DO*(O1(1:end-2)+O1(3:end)-2*O1(2:end-1))-GammaO(2:end-1)+KO*((Oex-O1(2:end-1))+abs(Oex-O1(2:end-1)))/2;
    Ab(1)=Ab(2); Ab(end)=Ab(end-1);    ub(:,1)=ub(:,2);ub(:,end)=ub(:,end-1); Ob(1)=Ob(2); Ob(end)=Ob(end-1);   
    %%
    u1=u+ub*dt/2;A1=A+Ab*dt/2;A1(A1<=0)=0;
    Gamma1=G0*sum(u1)./(1+(0.5./A1)).*((kaO./(1+((Oex-O1)/ktO).^5))+1-kaO);
    O1=O+Ob*dt/2;O1(O1<=0)=0;
    GammaO=GOx*sum(u1)./(1+10./O1);
        %f=ki'*A1;
        f=ki'*log((1+A1/Ki)./(1+A1/Ka));%
    % fO=kiO'*(log((1+O1/KiO1)./(1+O1/KaO1))-log((1+O1/KiO2)./(1+O1/KaO2)));
    g(:,2:end-1)=u1(:,2:end-1).*(f(:,3:end)-f(:,1:end-2))/2/dx;   g(:,1)=-g(:,2);g(:,end)=-g(:,end-1);
%      gO(:,2:end-1)=u1(:,2:end-1).*(fO(:,3:end)-fO(:,1:end-2))/2/dx;   gO(:,1)=-gO(:,2);gO(:,end)=-gO(:,end-1);
    uc(:,2:end-1)=Drho(:,2:end-1).*(u1(:,1:end-2)+u1(:,3:end)-2*u1(:,2:end-1))-(g(:,3:end)-g(:,1:end-2))/2/dx;%-(gO(:,3:end)-gO(:,1:end-2))/2/dx;
    Ac(2:end-1)=Dn*(A1(1:end-2)+A1(3:end)-2*A1(2:end-1))-Gamma1(2:end-1);
    Oc(2:end-1)=DO*(O1(1:end-2)+O1(3:end)-2*O1(2:end-1))-GammaO(2:end-1)+KO*((Oex-O1(2:end-1))+abs(Oex-O1(2:end-1)))/2;
     Ac(1)=Ac(2);Ac(end)=Ac(end-1);    uc(:,1)=uc(:,2);uc(:,end)=uc(:,end-1);Oc(1)=Oc(2); Oc(end)=Oc(end-1);   
    
    %%
    u1=u+uc*dt;A1=A+Ac*dt;A1(A1<=0)=0;
    Gamma1=G0*sum(u1)./(1+(0.5./A1)).*((kaO./(1+((Oex-O1)/ktO).^5))+1-kaO);
     O1=O+Oc*dt;O1(O1<=0)=0;
     GammaO=GOx*sum(u1)./(1+10./O1);
       %f=ki'*A1;
       f=ki'*log((1+A1/Ki)./(1+A1/Ka));%
    % fO=kiO'*(log((1+O1/KiO1)./(1+O1/KaO1))-log((1+O1/KiO2)./(1+O1/KaO2)));
    g(:,2:end-1)=u1(:,2:end-1).*(f(:,3:end)-f(:,1:end-2))/2/dx;   g(:,1)=-g(:,2);g(:,end)=-g(:,end-1);
%      gO(:,2:end-1)=u1(:,2:end-1).*(fO(:,3:end)-fO(:,1:end-2))/2/dx;   gO(:,1)=-gO(:,2);gO(:,end)=-gO(:,end-1);
    ud(:,2:end-1)=Drho(:,2:end-1).*(u1(:,1:end-2)+u1(:,3:end)-2*u1(:,2:end-1))-(g(:,3:end)-g(:,1:end-2))/2/dx;%-(gO(:,3:end)-gO(:,1:end-2))/2/dx;
    Ad(2:end-1)=Dn*(A1(1:end-2)+A1(3:end)-2*A1(2:end-1))-Gamma1(2:end-1);
    Od(2:end-1)=DO*(O1(1:end-2)+O1(3:end)-2*O1(2:end-1))-GammaO(2:end-1)+KO*((Oex-O1(2:end-1))+abs(Oex-O1(2:end-1)))/2;
    Ad(1)=Ad(2); Ad(end)=Ad(end-1);    ud(:,1)=ud(:,2);ud(:,end)=ud(:,end-1);Od(1)=Od(2); Od(end)=Od(end-1);   
    u=u+1/6*(ua+2*ub+2*uc+ud)*dt;
     A=A+1/6*(Aa+2*Ab+2*Ac+Ad)*dt;
     O=O+1/6*(Oa+2*Ob+2*Oc+Od)*dt;
     A(1)=A(2); A(end)=A(end-1);     u(:,1)=u(:,2);u(:,end)=u(:,end-1); 
     A(A<=0)=0;
     O(1)=O(2); O(end)=O(end-1); O(O<=0)=0;  
  
     indx=find(diff(sign(diff(sum(u))))==-2)+1;
     if ~isempty(indx)
        peakpos(t)=LX(indx(end));
     else 
         peakpos(t)=LX(1);
     end
       if (mod(t,120/dt)==0)
%          plot(LX,u1(:,1),'g');hold on;plot(LX,u1(:,2),'b');
%        for ii=1:length(CWBias)
%         plot(LX,u(ii,:),'color',CC(ii,:));hold on;

%        end%plot(n,'g');
        ut(:,:,tt)=u;
        At(:,tt)=A;
        Ot(:,tt)=O;
        tt=tt+1;
        subplot(2,1,1)
        plot(LX,sum(u(:,:))/max(sum(u(:,:))),'b');hold on;
        plot(LX,A/max(A),'r--'); plot(LX,O/max(O),'g--'); plot(LX,G0./(1+0.5./A).*(0.2+(O1/80).^2)./(1+(O1/80).^2)/G0,'c');
            subplot(2,1,2)
        plot(LX,sum(u(:,:)),'b');hold on;
         toc
     end
 
end
sum(sum(u(:,2:end-1)))
time=0:dt:(T-1)*dt;
p=polyfit(time(peakpos>2.5),peakpos(peakpos>2.5),1);
wavespeed=p(1)*1000;
final_profile.cell_density=u;
final_profile.asp=A;
final_profile.oxygen=O;
dyn_profile.cell_density=ut;
dyn_profile.Asp=At;
dyn_profile.O=Ot;

 end
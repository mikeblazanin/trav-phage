%% 
function [wavespeed_late final_profile dyn_profile dx]=chemo_oxygen_fun_v2_no_oxg_grad(Asp,Time)
%  dCWBias=0.005; CWBias=0.01:dCWBias:0.99;P = lognpdf(CWBias,log(0.295),0.308)/sum(lognpdf(CWBias,log(0.295),0.308))/dCWBias;
dCWBias=0.001; CWBias=0.001:dCWBias:0.99;P = lognpdf(CWBias,log(0.318),0.298)/sum(lognpdf(CWBias,log(0.318),0.298))/dCWBias;% define the original distribution; get from the experimental distribution, fit into log normal distribution. note the tails on both has slightly offset.  
dt=0.08;    %unit is s
dx=20;      % unit is um
L=13000/dx; %size of the channel 
LX=(0:dx:(round(L)*dx-dx))/1000; 
T=round(Time/dt);% simulation time
%% convert CWBias to chi and mu, according to our plos comp bio paper.
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
ki=theta*v^2*N*tau.*(lambdaRdiff./(1+lambdaR./lambdaT)./(lambdaR*theta+Drot))./(1+tau.*(lambdaR+Drot))/3;
%% preset parameters

Ki=0.3;      %uM
Ka=1000;    %uM
kaO=0.65; % asp consumption basel level
ktO=60; % asp oxy dependent threshold
cell_density=zeros(length(CWBias),round(L)); %cell density 
Drho=Diff'/dx^2*ones(1,round(L)); % diffusion constant for every 
Asp_field=ones(1,round(L))*Asp; %Asp_field in uM
Oex=250;% external oxy level in uM
Oxy_field=ones(1,round(L))*Oex;% oxy field
KO=0.02; %%PDMS oxy transfer rate
GOx=8*4/60;% oxy consumption rate
Dn=600/dx^2;%Asp method for the determination of diffusion coefficients for small molecules in aqueous solution, Analytical Biochemistry,Volume 166, Issue 2, 1 November 1987, Pages 335-341
DO=2500/dx^2;%http://compost.css.cornell.edu/oxygen/oxygen.diff.water.html
G0=8.5*0.84/60;% 10uM/min for OD 1 cells.


% define temp fields for 4th-RK 
cell_density_A=cell_density;
cell_density_B=cell_density;
cell_density_C=cell_density;
cell_density_D=cell_density; 
Asp_field_A=Asp_field;
Asp_field_B=Asp_field;
Asp_field_C=Asp_field;
Asp_field_D=Asp_field;
Oxy_field_A=Oxy_field;
Oxy_field_B=Oxy_field;
Oxy_field_C=Oxy_field;
Oxy_field_D=Oxy_field;
tt=1;
g=cell_density;

%% initial conditions
DL=1000/dx; V=1000/dx;
for x=1:round(DL)
    ux(x)=1;%exp(-(x)^2/V^2);
end
ux=ux./sum(ux);
cell_density(:,2:round(DL)+1)=P'*ux*500*50/dx*dCWBias; %unit is OD1,~6*8.4*1.4*dx cells in 600um*14um*dx %total cell: 0.7*10^9*(600*14*10^4+900*14*900+?)*10^-12=15~20*10^4 cells
cell_density(:,1)=cell_density(:,2);
Asp_field=Asp*(1-1./(1+(LX*1000/(V*dx*3)).^2.5));% initial the Asp field
Oxy_field=Oex*(1-1./(1+(LX*1000/(V*dx*6)).^3)); % initial the oxy field



%%
for t=1:T
      %% 4th order RK method
    cell_density_temp=cell_density;
    Asp_temp=Asp_field;
    Asp_temp(Asp_temp<=0)=0;
    Oxy_field_temp=Oxy_field; 
    Oxy_field_temp(Oxy_field_temp<=0)=0;
    Gamma1=G0*sum(cell_density_temp)./(1+(0.5./Asp_temp)).*((kaO./(1+((Oex-Oxy_field_temp)/ktO).^5))+1-kaO);% Asp consumption
    GammaO=GOx*sum(cell_density_temp)./(1+10./Oxy_field_temp); %oxy consumption
    f=ki'*log((1+Asp_temp/Ki)./(1+Asp_temp/Ka));% 
    g(:,2:end-1)=cell_density_temp(:,2:end-1).*(f(:,3:end)-f(:,1:end-2))/2/dx;   
    g(:,1)=-g(:,2);
    g(:,end)=-g(:,end-1);
    cell_density_A(:,2:end-1)=Drho(:,2:end-1).*(cell_density_temp(:,1:end-2)+cell_density_temp(:,3:end)-2*cell_density_temp(:,2:end-1))-(g(:,3:end)-g(:,1:end-2))/2/dx;
    Asp_field_A(2:end-1)=Dn*(Asp_temp(1:end-2)+Asp_temp(3:end)-2*Asp_temp(2:end-1))-Gamma1(2:end-1);
    Oxy_field_A(2:end-1)=DO*(Oxy_field_temp(1:end-2)+Oxy_field_temp(3:end)-2*Oxy_field_temp(2:end-1))-GammaO(2:end-1)+KO*((Oex-Oxy_field_temp(2:end-1))+abs(Oex-Oxy_field_temp(2:end-1)))/2;
         % Boundary conditions
    Asp_field_A(1)=Asp_field_A(2); 
    Asp_field_A(end)=Asp_field_A(end-1);    
    cell_density_A(:,1)=cell_density_A(:,2);
    cell_density_A(:,end)=cell_density_A(:,end-1); 
    Oxy_field_A(1)=Oxy_field_A(2); 
    Oxy_field_A(end)=Oxy_field_A(end-1);   
 
    %
    cell_density_temp=cell_density+cell_density_A*dt/2;
    Asp_temp=Asp_field+Asp_field_A*dt/2;Asp_temp(Asp_temp<=0)=0;
    Gamma1=G0*sum(cell_density_temp)./(1+(0.5./Asp_temp)).*((kaO./(1+((Oex-Oxy_field_temp)/ktO).^5))+1-kaO);
    Oxy_field_temp=Oxy_field+Oxy_field_A*dt/2;Oxy_field_temp(Oxy_field_temp<=0)=0;
    GammaO=GOx*sum(cell_density_temp)./(1+10./Oxy_field_temp);
    f=ki'*log((1+Asp_temp/Ki)./(1+Asp_temp/Ka));%
    g(:,2:end-1)=cell_density_temp(:,2:end-1).*(f(:,3:end)-f(:,1:end-2))/2/dx;   
    g(:,1)=-g(:,2);
    g(:,end)=-g(:,end-1);
    cell_density_B(:,2:end-1)=Drho(:,2:end-1).*(cell_density_temp(:,1:end-2)+cell_density_temp(:,3:end)-2*cell_density_temp(:,2:end-1))-(g(:,3:end)-g(:,1:end-2))/2/dx;
    Asp_field_B(2:end-1)=Dn*(Asp_temp(1:end-2)+Asp_temp(3:end)-2*Asp_temp(2:end-1))-Gamma1(2:end-1);
    Oxy_field_B(2:end-1)=DO*(Oxy_field_temp(1:end-2)+Oxy_field_temp(3:end)-2*Oxy_field_temp(2:end-1))-GammaO(2:end-1)+KO*((Oex-Oxy_field_temp(2:end-1))+abs(Oex-Oxy_field_temp(2:end-1)))/2;
         % Boundary conditions
    Asp_field_B(1)=Asp_field_B(2); 
    Asp_field_B(end)=Asp_field_B(end-1);    
    cell_density_B(:,1)=cell_density_B(:,2);
    cell_density_B(:,end)=cell_density_B(:,end-1); 
    Oxy_field_B(1)=Oxy_field_B(2); 
    Oxy_field_B(end)=Oxy_field_B(end-1);   
    %
    cell_density_temp=cell_density+cell_density_B*dt/2;
    Asp_temp=Asp_field+Asp_field_B*dt/2;
    Asp_temp(Asp_temp<=0)=0;
    Gamma1=G0*sum(cell_density_temp)./(1+(0.5./Asp_temp)).*((kaO./(1+((Oex-Oxy_field_temp)/ktO).^5))+1-kaO);
    Oxy_field_temp=Oxy_field+Oxy_field_B*dt/2;
    Oxy_field_temp(Oxy_field_temp<=0)=0;
    GammaO=GOx*sum(cell_density_temp)./(1+10./Oxy_field_temp);
    f=ki'*log((1+Asp_temp/Ki)./(1+Asp_temp/Ka));%
    g(:,2:end-1)=cell_density_temp(:,2:end-1).*(f(:,3:end)-f(:,1:end-2))/2/dx;   
    g(:,1)=-g(:,2);
    g(:,end)=-g(:,end-1);
    cell_density_C(:,2:end-1)=Drho(:,2:end-1).*(cell_density_temp(:,1:end-2)+cell_density_temp(:,3:end)-2*cell_density_temp(:,2:end-1))-(g(:,3:end)-g(:,1:end-2))/2/dx;
    Asp_field_C(2:end-1)=Dn*(Asp_temp(1:end-2)+Asp_temp(3:end)-2*Asp_temp(2:end-1))-Gamma1(2:end-1);
    Oxy_field_C(2:end-1)=DO*(Oxy_field_temp(1:end-2)+Oxy_field_temp(3:end)-2*Oxy_field_temp(2:end-1))-GammaO(2:end-1)+KO*((Oex-Oxy_field_temp(2:end-1))+abs(Oex-Oxy_field_temp(2:end-1)))/2;
         % Boundary conditions
    Asp_field_C(1)=Asp_field_C(2);
    Asp_field_C(end)=Asp_field_C(end-1);    
    cell_density_C(:,1)=cell_density_C(:,2);
    cell_density_C(:,end)=cell_density_C(:,end-1);
    Oxy_field_C(1)=Oxy_field_C(2); 
    Oxy_field_C(end)=Oxy_field_C(end-1);   
    
    %%
    cell_density_temp=cell_density+cell_density_D*dt;Asp_temp=Asp_field+Asp_field_D*dt;Asp_temp(Asp_temp<=0)=0;
    Gamma1=G0*sum(cell_density_temp)./(1+(0.5./Asp_temp)).*((kaO./(1+((Oex-Oxy_field_temp)/ktO).^5))+1-kaO);
     Oxy_field_temp=Oxy_field+Oxy_field_D*dt;Oxy_field_temp(Oxy_field_temp<=0)=0;
     GammaO=GOx*sum(cell_density_temp)./(1+10./Oxy_field_temp);
     f=ki'*log((1+Asp_temp/Ki)./(1+Asp_temp/Ka));%
    g(:,2:end-1)=cell_density_temp(:,2:end-1).*(f(:,3:end)-f(:,1:end-2))/2/dx;   
    g(:,1)=-g(:,2);
    g(:,end)=-g(:,end-1);
    cell_density_D(:,2:end-1)=Drho(:,2:end-1).*(cell_density_temp(:,1:end-2)+cell_density_temp(:,3:end)-2*cell_density_temp(:,2:end-1))-(g(:,3:end)-g(:,1:end-2))/2/dx;
    Asp_field_D(2:end-1)=Dn*(Asp_temp(1:end-2)+Asp_temp(3:end)-2*Asp_temp(2:end-1))-Gamma1(2:end-1);
    Oxy_field_D(2:end-1)=DO*(Oxy_field_temp(1:end-2)+Oxy_field_temp(3:end)-2*Oxy_field_temp(2:end-1))-GammaO(2:end-1)+KO*((Oex-Oxy_field_temp(2:end-1))+abs(Oex-Oxy_field_temp(2:end-1)))/2;
     % Boundary conditions
    Asp_field_D(1)=Asp_field_D(2); 
    Asp_field_D(end)=Asp_field_D(end-1);    
    cell_density_D(:,1)=cell_density_D(:,2);
    cell_density_D(:,end)=cell_density_D(:,end-1);
    Oxy_field_D(1)=Oxy_field_D(2); 
    Oxy_field_D(end)=Oxy_field_D(end-1);  
    %% update the the cell density, asp, oxy fields
     cell_density=cell_density+1/6*(cell_density_A+2*cell_density_B+2*cell_density_C+cell_density_D)*dt;
     Asp_field=Asp_field+1/6*(Asp_field_A+2*Asp_field_B+2*Asp_field_C+Asp_field_D)*dt;
     Oxy_field=Oxy_field+1/6*(Oxy_field_A+2*Oxy_field_B+2*Oxy_field_C+Oxy_field_D)*dt;
     % Boundary conditions
     Asp_field(1)=Asp_field(2); 
     Asp_field(end)=Asp_field(end-1);     
     Asp_field(Asp_field<=0)=0;
     cell_density(:,1)=cell_density(:,2);
     cell_density(:,end)=cell_density(:,end-1); 
     Oxy_field(1)=Oxy_field(2); 
     Oxy_field(end)=Oxy_field(end-1); 
     Oxy_field(Oxy_field<=0)=0;  
  %% do plots
     indx=find(diff(sign(diff(sum(cell_density))))==-2)+1;
     if ~isempty(indx)
        peakpos(t)=LX(indx(end));
     else 
         peakpos(t)=LX(1);
     end
%       if t>10*60/dt
       if (mod(t,60/dt)==0)
%          plot(LX,cell_density_temp(:,1),'g');hold on;plot(LX,cell_density_temp(:,2),'b');
%        for ii=1:length(CWBias)
%         plot(LX,cell_density(ii,:),'color',CC(ii,:));hold on;

%        end%plot(n,'g');
        cell_density_out(:,:,tt)=cell_density;
        Asp_field_out(:,tt)=Asp_field;
        Oxy_field_out(:,tt)=Oxy_field;
        tt=tt+1;
        subplot(2,1,1)
        plot(LX,sum(cell_density(:,:))/max(sum(cell_density(:,:))),'b');hold on;
        plot(LX,Asp_field/max(Asp_field),'r--'); plot(LX,Oxy_field/max(Oxy_field),'g--'); plot(LX,G0./(1+0.5./Asp_field).*(0.2+(Oxy_field_temp/80).^2)./(1+(Oxy_field_temp/80).^2)/G0,'c');
            subplot(2,1,2)
        plot(LX,sum(cell_density(:,:)),'b');hold on;
        
          toc
%      end
     end
end
sum(sum(cell_density(:,2:end-1)));
time=0:dt:(T-1)*dt;
p=polyfit(time(peakpos>6),peakpos(peakpos>6),1);

wavespeed_late=p(1)*1000;

final_profile.cell_density=cell_density;
final_profile.asp=Asp_field;
final_profile.oxygen=Oxy_field;
dyn_profile.cell_density=cell_density_out;
dyn_profile.Asp_field=Asp_field_out;
dyn_profile.Oxy_field=Oxy_field_out;

 end
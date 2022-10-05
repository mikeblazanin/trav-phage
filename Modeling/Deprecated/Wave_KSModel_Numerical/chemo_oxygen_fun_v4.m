%% calculate the chemotactic coefficients based on the PLOS Comp. 
 function [final_profile wavespeed dyn_profile dx]=chemo_oxygen_fun_v4(Asp,Time,Time0)
%  dCWBias=0.005; CWBias=0.01:dCWBias:0.99;P = lognpdf(CWBias,log(0.295),0.308)/sum(lognpdf(CWBias,log(0.295),0.308))/dCWBias;
CWBias_exp=[0,0.0100000000000000,0.0200000000000000,0.0300000000000000,0.0400000000000000,0.0500000000000000,0.0600000000000000,0.0700000000000000,0.0800000000000000,0.0900000000000000,0.100000000000000,0.110000000000000,0.120000000000000,0.130000000000000,0.140000000000000,0.150000000000000,0.160000000000000,0.170000000000000,0.180000000000000,0.190000000000000,0.200000000000000,0.210000000000000,0.220000000000000,0.230000000000000,0.240000000000000,0.250000000000000,0.260000000000000,0.270000000000000,0.280000000000000,0.290000000000000,0.300000000000000,0.310000000000000,0.320000000000000,0.330000000000000,0.340000000000000,0.350000000000000,0.360000000000000,0.370000000000000,0.380000000000000,0.390000000000000,0.400000000000000,0.410000000000000,0.420000000000000,0.430000000000000,0.440000000000000,0.450000000000000,0.460000000000000,0.470000000000000,0.480000000000000,0.490000000000000,0.500000000000000,0.510000000000000,0.520000000000000,0.530000000000000,0.540000000000000,0.550000000000000,0.560000000000000,0.570000000000000,0.580000000000000,0.590000000000000,0.600000000000000,0.610000000000000,0.620000000000000,0.630000000000000,0.640000000000000,0.650000000000000,0.660000000000000,0.670000000000000,0.680000000000000,0.690000000000000,0.700000000000000];
P_exp=[0,0.00302823490236943,0.00498640060510877,0.00740654246062159,0.0129576724735464,0.0282364445536801,0.0342746634745129,0.0605832532805051,0.0797351468975508,0.105863629784042,0.144957679208212,0.222281447976504,0.281608928805501,0.373987712933409,0.512065751136908,0.691152325282671,0.919556834551468,1.19485305314468,1.51155357679430,1.88448868628358,2.26163802431774,2.67050119556754,3.04502589236555,3.42156794182481,3.78827243915744,4.11858437055020,4.26770099008759,4.46704480752076,4.50323811436766,4.50544966346908,4.42742856648661,4.31612752213579,4.11316467718262,3.94529362322380,3.65109481682462,3.23848249620281,2.98678093597204,2.69627429599337,2.36626462681502,2.10807472740901,1.95992530735338,1.71658120173055,1.53144380527692,1.43826939988074,1.27455605283792,1.14836389445003,1.01135546697596,0.907225653941910,0.780385038989034,0.685196631521439,0.645362573778435,0.620271119730271,0.579787157033761,0.525204939995039,0.523233085565830,0.466362640857505,0.445131198472857,0.399401650650870,0.399447904195710,0.371516048640824,0.380986255044554,0.393076006538824,0.375160618566594,0.357271202286408,0.341376560968248,0.322917046544370,0.291493170711647,0.311004430859554,0.290374229335908,0.272975724047314,0.162176873861474];
dCWBias=0.005;CWBias=0.001:dCWBias:CWBias_exp(end);
P=interp1(CWBias_exp,P_exp,CWBias,'spline');P=P./sum(P)./dCWBias;
%      CWBias=[0.2 0.2];P=[0.5 0.5];dCWBias=1;
% CC=winter(length(CWBias));
dt=0.5;    %s
dx=50;      %um

%%
%     dt=0.5;    %s
% dx=50;      %um
g1=28;g2=g1;
Kd      =  2.9;
w       =  3.8;
Yp=g2*Kd./(g2-g1/2+log(1./CWBias-1))-Kd;
deltaG  =  g1/4-g2/2*(Yp./(Kd+Yp));
lambdaT   =  w*exp(deltaG);
lambdaR  =  w*exp(-deltaG);
CWBias  =  ((lambdaR)./((lambdaT)+(lambdaR)));
 alpha=5.1138;

Drot=0.062*2;
v=28*4/pi;%2D run speed 26um/s
 theta=1-0.1564;
Diff=v^2*(1-CWBias)./3./(lambdaR*theta+Drot);
ki=22*Diff./(1+0.05./CWBias);%.*(0.995./CWBias)./(1+0.995./CWBias);
%%
% cwbin=[0,0.0250000000000000,0.0750000000000000,0.125000000000000,0.175000000000000,0.225000000000000,0.275000000000000,0.325000000000000,0.375000000000000,0.425000000000000,0.475000000000000,0.525000000000000,0.575000000000000,0.7];
% ki_exp=[0 7902.09549303570,8073.05214050945,5981.70231869196,4355.52930404505,3103.66077946146,2174.74772106217,1557.66888958829,1020.05508365524,900,766.694837600473,419.861936116645,335,150];
% % % dCWBias=0.001;CWBias=CWBias_exp(1):dCWBias:CWBias_exp(end);
%  ki=interp1(cwbin,ki_exp,CWBias,'spline');
% cwbin=[0	0.01	0.03	0.05	0.07	0.09	0.11	0.13	0.15	0.17	0.19	0.21	0.23	0.25	0.27	0.29	0.31	0.33	0.35	0.37	0.39	0.41	0.43	0.45	0.47	0.49	0.51	0.53	0.55	0.57	0.59	0.61	0.63	0.65	0.67	0.69	0.725	0.775	0.825	0.875	0.925	0.975];
% mu_exp=[1800	1721.946291	945.6181097	698.2559545	575.7232203	500.0422885	437.2564495	363.0177159	320.3600442	275.2751298	244.1193479	210.2766391	187.8269693	166.9712219	145.0543619	130.1950861	116.6899906	103.8049972	92.39097476	83.59047223	76.46145915	66.06179514	63.52649191	54.78469659	50.16759063	46.96056098	41.5261955	36.43790868	32.42582081	32.80620252	30.65181085	23.71780178	25.06033671	22.31129005	18.52295702	16.05250083	14.71569087	11.39640718	9.6552	7.1062	3.937230494	2];
% % cwbin=[0	0.025	0.075	0.125	0.175	0.225	0.275	0.325	0.375	0.425	0.475	0.525	0.575	0.625	0.675	0.725	0.775	0.825	0.875	0.925	0.975];
% % mu_exp=[1100.377077	1060.377077	917.939342	627.808247	430.5267209	292.8775278	205.0181211	146.1547262	101.6368463	72.43995219	53.18062237	43.84310308	35.5219	25.7262	18.861936	14.71569087	11.39640718	9.6552	7.1062	3.937230494	2];
% % dCWBias=0.001;CWBias=CWBias_exp(1):dCWBias:CWBias_exp(end);
% Diff=interp1(cwbin,mu_exp,CWBias,'spline');
% ki=30*Diff./(1+0.12./CWBias);%.*(0.995./CWBias)./(1+0.995./CWBias);
% % % kiO=ki/5000;
%%   

Ki=3.5;      %uM
Ka=1000;    %uM
% KiO1=40;% KiO2=150;% KaO1=330;% KaO2=10000;% Studies of bacterial aerotaxis in a microfluidic device 15% transition from attractant to reppellent%  Asp=100;   %initial Asp conc in units of uM
% kaO=0.25;
% ktO=k;
% kappa=1.5;
kappa=1.;
alpha0=0.27;
kA=0.5;
L=20000/dx;
LX=(0:dx:(L*dx-dx))/1000;
u=zeros(length(CWBias),round(L));
Drho=Diff'/dx^2*ones(1,round(L));
A=ones(1,round(L))*Asp;
Oex=250;
O=ones(1,round(L))*Oex;
% norm_alpha_O=(1-kaO./(1+(Oex./ktO).^kappa));
KO=0.02;
GOx=7*8.4/60;
  T=round(Time/dt);%
  T1=round(Time0/dt);
Dn=500/dx^2;%A method for the determination of diffusion coefficients for small molecules in aqueous solution, Analytical Biochemistry,Volume 166, Issue 2, 1 November 1987, Pages 335-341
DO=2500/dx^2;%http://compost.css.cornell.edu/oxygen/oxygen.diff.water.html
G0=9.3*0.84/60;% 9uM/min for OD 1 cells.
DL=800/dx; %V=450/dx;%
V=500/dx;
for x=1:round(2*DL)
    ux(x)=exp(-(x)^2/(1.*DL)^2);
end
ux=ux./sum(ux);
u(:,2:round(2*DL)+1)=P'*ux*500*47.5/dx*dCWBias; %unit is OD1,~6*8.4*1.4*dx cells in 600um*14um*dx %total cell: 0.7*10^9*(600*14*10^4+900*14*900+?)*10^-12=15~20*10^4 cells
u(:,1)=u(:,2);
A=Asp*(1-1./(1+(LX*1000/(V*dx*3)).^3));
% A=Asp*(1-1./(1+(LX*1000/(V*dx*2)).^0));
% O=Oex*(1-0.8./(1+(LX*1000/(V*dx*2)).^6));
ua=u;ub=u;uc=u;ud=u;Aa=A;Ab=A;Ac=A;Ad=A;Oa=O;Ob=O;Oc=O;Od=O;


%%
tt=1;
g=u;grad=A;gO=u;
%%
for t=1:T1
    u1=u;A1=A;A1(A1<=0)=0;
    O1=O; O1(O1<=0)=0;
    O1_sense=O1;%O1_sense=smooth(O1,span)';
    Gamma1=G0*sum(u1)./(1+(kA./A1)).*((1-alpha0)*(O1_sense./Oex).^kappa+alpha0);%(kaO./(1+((Oex-O)/ktO).^5))+1-kaO;
    GammaO=GOx*sum(u1)./(1+1./O1);


    %f=ki'*A1;
    Aa(2:end-1)=Dn*(A1(1:end-2)+A1(3:end)-2*A1(2:end-1))-Gamma1(2:end-1);
    Oa(2:end-1)=DO*(O1(1:end-2)+O1(3:end)-2*O1(2:end-1))-GammaO(2:end-1)+KO*((Oex-O1_sense(2:end-1))+abs(Oex-O1_sense(2:end-1)))/2;
    Aa(1)=Aa(2); Aa(end)=Aa(end-1);    Oa(1)=Oa(2); Oa(end)=Oa(end-1);   
 
    %%
    A1=A+Aa*dt/2;A1(A1<=0)=0;
    O1=O+Oa*dt/2;O1(O1<=0)=0;
    O1_sense=O1;%O1_sense=smooth(O1,span)';
    Gamma1=G0*sum(u1)./(1+(kA./A1)).*((1-alpha0)*(O1_sense./Oex).^kappa+alpha0);
    GammaO=GOx*sum(u1)./(1+1./O1);
    %f=ki'*A1;
    Ab(2:end-1)=Dn*(A1(1:end-2)+A1(3:end)-2*A1(2:end-1))-Gamma1(2:end-1);
    Ob(2:end-1)=DO*(O1(1:end-2)+O1(3:end)-2*O1(2:end-1))-GammaO(2:end-1)+KO*((Oex-O1_sense(2:end-1))+abs(Oex-O1_sense(2:end-1)))/2;
    Ab(1)=Ab(2); Ab(end)=Ab(end-1);    Ob(1)=Ob(2); Ob(end)=Ob(end-1);   
    %%
    A1=A+Ab*dt/2;A1(A1<=0)=0;
    O1=O+Ob*dt/2;O1(O1<=0)=0;
    O1_sense=O1;%O1_sense=smooth(O1,span)';
    Gamma1=G0*sum(u1)./(1+(kA./A1)).*((1-alpha0)*(O1_sense./Oex).^kappa+alpha0);
    GammaO=GOx*sum(u1)./(1+1./O1);
    Ac(2:end-1)=Dn*(A1(1:end-2)+A1(3:end)-2*A1(2:end-1))-Gamma1(2:end-1);
    Oc(2:end-1)=DO*(O1(1:end-2)+O1(3:end)-2*O1(2:end-1))-GammaO(2:end-1)+KO*((Oex-O1_sense(2:end-1))+abs(Oex-O1_sense(2:end-1)))/2;
     Ac(1)=Ac(2);Ac(end)=Ac(end-1);    Oc(1)=Oc(2); Oc(end)=Oc(end-1);   
    
    %%
    A1=A+Ac*dt;A1(A1<=0)=0;
    O1=O+Oc*dt;O1(O1<=0)=0; 
    O1_sense=O1;%O1_sense=smooth(O1,span)';
    Gamma1=G0*sum(u1)./(1+(kA./A1)).*((1-alpha0)*(O1_sense./Oex).^kappa+alpha0);
    GammaO=GOx*sum(u1)./(1+1./O1);
    Ad(2:end-1)=Dn*(A1(1:end-2)+A1(3:end)-2*A1(2:end-1))-Gamma1(2:end-1);
    Od(2:end-1)=DO*(O1(1:end-2)+O1(3:end)-2*O1(2:end-1))-GammaO(2:end-1)+KO*((Oex-O1_sense(2:end-1))+abs(Oex-O1_sense(2:end-1)))/2;
    Ad(1)=Ad(2); Ad(end)=Ad(end-1);   Od(1)=Od(2); Od(end)=Od(end-1);   
    A=A+1/6*(Aa+2*Ab+2*Ac+Ad)*dt;
    O=O+1/6*(Oa+2*Ob+2*Oc+Od)*dt;
    A(1)=A(2); A(end)=A(end-1);   
    A(A<=0)=0;
    O(1)=O(2); O(end)=O(end-1); O(O<=0)=0;  
     
             
end
%%
for t=1:T
    u1=u;A1=A;A1(A1<=0)=0;
    O1=O; O1(O1<=0)=0;
    O1_sense=O1;%O1_sense=smooth(O1,span)';
    Gamma1=G0*sum(u1)./(1+(kA./A1)).*((1-alpha0)*(O1_sense./Oex).^kappa+alpha0);%(kaO./(1+((Oex-O)/ktO).^5))+1-kaO;
    GammaO=GOx*sum(u1)./(1+1./O1);


    %f=ki'*A1;
     f=ki'*log((1+A1/Ki)./(1+A1/Ka));%
    % fO=kiO'*(log((1+O1/KiO1)./(1+O1/KaO1))-log((1+O1/KiO2)./(1+O1/KaO2)));
    g(:,2:end-1)=u1(:,2:end-1).*(f(:,3:end)-f(:,1:end-2))/2/dx;   g(:,1)=-g(:,2);g(:,end)=-g(:,end-1);
%      gO(:,2:end-1)=u1(:,2:end-1).*(fO(:,3:end)-fO(:,1:end-2))/2/dx;   gO(:,1)=-gO(:,2);gO(:,end)=-gO(:,end-1);
    ua(:,2:end-1)=Drho(:,2:end-1).*(u1(:,1:end-2)+u1(:,3:end)-2*u1(:,2:end-1))-(g(:,3:end)-g(:,1:end-2))/2/dx;%-(gO(:,3:end)-gO(:,1:end-2))/2/dx;
    Aa(2:end-1)=Dn*(A1(1:end-2)+A1(3:end)-2*A1(2:end-1))-Gamma1(2:end-1);
    Oa(2:end-1)=DO*(O1(1:end-2)+O1(3:end)-2*O1(2:end-1))-GammaO(2:end-1)+KO*((Oex-O1_sense(2:end-1))+abs(Oex-O1_sense(2:end-1)))/2;
    Aa(1)=Aa(2); Aa(end)=Aa(end-1);    ua(:,1)=ua(:,2);ua(:,end)=ua(:,end-1);  Oa(1)=Oa(2); Oa(end)=Oa(end-1);   
 
    %%
    u1=u+ua*dt/2;A1=A+Aa*dt/2;A1(A1<=0)=0;
    
    O1=O+Oa*dt/2;O1(O1<=0)=0;
     O1_sense=O1;%O1_sense=smooth(O1,span)';
    Gamma1=G0*sum(u1)./(1+(kA./A1)).*((1-alpha0)*(O1_sense./Oex).^kappa+alpha0);%(kaO./(1+((Oex-O)/ktO).^5))+1-kaO;
    GammaO=GOx*sum(u1)./(1+1./O1);
    %f=ki'*A1;
    f=ki'*log((1+A1/Ki)./(1+A1/Ka));%
    % fO=kiO'*(log((1+O1/KiO1)./(1+O1/KaO1))-log((1+O1/KiO2)./(1+O1/KaO2)));
    g(:,2:end-1)=u1(:,2:end-1).*(f(:,3:end)-f(:,1:end-2))/2/dx;   g(:,1)=-g(:,2);g(:,end)=-g(:,end-1);
%      gO(:,2:end-1)=u1(:,2:end-1).*(fO(:,3:end)-fO(:,1:end-2))/2/dx;   gO(:,1)=-gO(:,2);gO(:,end)=-gO(:,end-1);
    ub(:,2:end-1)=Drho(:,2:end-1).*(u1(:,1:end-2)+u1(:,3:end)-2*u1(:,2:end-1))-(g(:,3:end)-g(:,1:end-2))/2/dx;%-(gO(:,3:end)-gO(:,1:end-2))/2/dx;
    Ab(2:end-1)=Dn*(A1(1:end-2)+A1(3:end)-2*A1(2:end-1))-Gamma1(2:end-1);
    Ob(2:end-1)=DO*(O1(1:end-2)+O1(3:end)-2*O1(2:end-1))-GammaO(2:end-1)+KO*((Oex-O1_sense(2:end-1))+abs(Oex-O1_sense(2:end-1)))/2;
    Ab(1)=Ab(2); Ab(end)=Ab(end-1);    ub(:,1)=ub(:,2);ub(:,end)=ub(:,end-1); Ob(1)=Ob(2); Ob(end)=Ob(end-1);   
    %%
    u1=u+ub*dt/2;A1=A+Ab*dt/2;A1(A1<=0)=0;
    
    O1=O+Ob*dt/2;O1(O1<=0)=0; O1_sense=O1;%O1_sense=smooth(O1,span)';
    Gamma1=G0*sum(u1)./(1+(kA./A1)).*((1-alpha0)*(O1_sense./Oex).^kappa+alpha0);%(kaO./(1+((Oex-O)/ktO).^5))+1-kaO;
    GammaO=GOx*sum(u1)./(1+1./O1);
        f=ki'*log((1+A1/Ki)./(1+A1/Ka));%
    % fO=kiO'*(log((1+O1/KiO1)./(1+O1/KaO1))-log((1+O1/KiO2)./(1+O1/KaO2)));
    g(:,2:end-1)=u1(:,2:end-1).*(f(:,3:end)-f(:,1:end-2))/2/dx;   g(:,1)=-g(:,2);g(:,end)=-g(:,end-1);
%      gO(:,2:end-1)=u1(:,2:end-1).*(fO(:,3:end)-fO(:,1:end-2))/2/dx;   gO(:,1)=-gO(:,2);gO(:,end)=-gO(:,end-1);
    uc(:,2:end-1)=Drho(:,2:end-1).*(u1(:,1:end-2)+u1(:,3:end)-2*u1(:,2:end-1))-(g(:,3:end)-g(:,1:end-2))/2/dx;%-(gO(:,3:end)-gO(:,1:end-2))/2/dx;
    Ac(2:end-1)=Dn*(A1(1:end-2)+A1(3:end)-2*A1(2:end-1))-Gamma1(2:end-1);
    Oc(2:end-1)=DO*(O1(1:end-2)+O1(3:end)-2*O1(2:end-1))-GammaO(2:end-1)+KO*((Oex-O1_sense(2:end-1))+abs(Oex-O1_sense(2:end-1)))/2;
     Ac(1)=Ac(2);Ac(end)=Ac(end-1);    uc(:,1)=uc(:,2);uc(:,end)=uc(:,end-1);Oc(1)=Oc(2); Oc(end)=Oc(end-1);   
    
    %%
    u1=u+uc*dt;A1=A+Ac*dt;A1(A1<=0)=0;
    
     O1=O+Oc*dt;O1(O1<=0)=0; 
     O1_sense=O1;%O1_sense=smooth(O1,span)';
    Gamma1=G0*sum(u1)./(1+(kA./A1)).*((1-alpha0)*(O1_sense./Oex).^kappa+alpha0);%(kaO./(1+((Oex-O)/ktO).^5))+1-kaO;
    GammaO=GOx*sum(u1)./(1+1./O1);
       %f=ki'*A1;
       f=ki'*log((1+A1/Ki)./(1+A1/Ka));%
    % fO=kiO'*(log((1+O1/KiO1)./(1+O1/KaO1))-log((1+O1/KiO2)./(1+O1/KaO2)));
    g(:,2:end-1)=u1(:,2:end-1).*(f(:,3:end)-f(:,1:end-2))/2/dx;   g(:,1)=-g(:,2);g(:,end)=-g(:,end-1);
%      gO(:,2:end-1)=u1(:,2:end-1).*(fO(:,3:end)-fO(:,1:end-2))/2/dx;   gO(:,1)=-gO(:,2);gO(:,end)=-gO(:,end-1);
    ud(:,2:end-1)=Drho(:,2:end-1).*(u1(:,1:end-2)+u1(:,3:end)-2*u1(:,2:end-1))-(g(:,3:end)-g(:,1:end-2))/2/dx;%-(gO(:,3:end)-gO(:,1:end-2))/2/dx;
    Ad(2:end-1)=Dn*(A1(1:end-2)+A1(3:end)-2*A1(2:end-1))-Gamma1(2:end-1);
    Od(2:end-1)=DO*(O1(1:end-2)+O1(3:end)-2*O1(2:end-1))-GammaO(2:end-1)+KO*((Oex-O1_sense(2:end-1))+abs(Oex-O1_sense(2:end-1)))/2;
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
       if (mod(t,60/dt)==0)
%          plot(LX,u1(:,1),'g');hold on;plot(LX,u1(:,2),'b');
%         for ii=1:length(CWBias)
%         plot(LX,u(ii,:),'color',CC(ii,:));hold on;

%        end%plot(n,'g');
        ut(:,:,tt)=u;
        At(:,tt)=A;
        Ot(:,tt)=O;
        tt=tt+1;
%         subplot(2,1,1)
%         plot(LX,sum(u(:,:))/max(sum(u(:,:))),'b');hold on;
%         plot(LX,A/max(A),'r--'); plot(LX,O/max(O),'g--'); plot(LX,G0./(1+0.5./A).*(kaO+(O/ktO).^kappa)./(1+(O/ktO).^kappa)/G0,'c');
%             subplot(2,1,2)
%         plot(LX,sum(u(:,:)),'b');hold on;
%          toc
     end
 
end
sum(sum(u(:,2:end-1)))
time=0:dt:(T-1)*dt;
p=polyfit(time(peakpos>2.5),peakpos(peakpos>2.5),1);
wavespeed=p(1)*1000/1000*60;
final_profile.cell_density=u;
final_profile.asp=A;
final_profile.oxygen=O;
final_profile.peakpos=peakpos;
dyn_profile.cell_density=ut;
dyn_profile.Asp=At;
dyn_profile.O=Ot;

 end
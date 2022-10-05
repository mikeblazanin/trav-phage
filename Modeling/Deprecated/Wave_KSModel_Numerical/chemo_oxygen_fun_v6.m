%% calculate the chemotactic coefficients based on the PLOS Comp. 
 function [final_profile wavespeed dyn_profile dx]=chemo_oxygen_fun_v6(Asp,Time,Time0)
%  dCWBias=0.005; CWBias=0.01:dCWBias:0.99;P = lognpdf(CWBias,log(0.295),0.308)/sum(lognpdf(CWBias,log(0.295),0.308))/dCWBias;
CWBias_exp=[0,0.0100000000000000,0.0200000000000000,0.0300000000000000,0.0400000000000000,0.0500000000000000,0.0600000000000000,0.0700000000000000,0.0800000000000000,0.0900000000000000,0.100000000000000,0.110000000000000,0.120000000000000,0.130000000000000,0.140000000000000,0.150000000000000,0.160000000000000,0.170000000000000,0.180000000000000,0.190000000000000,0.200000000000000,0.210000000000000,0.220000000000000,0.230000000000000,0.240000000000000,0.250000000000000,0.260000000000000,0.270000000000000,0.280000000000000,0.290000000000000,0.300000000000000,0.310000000000000,0.320000000000000,0.330000000000000,0.340000000000000,0.350000000000000,0.360000000000000,0.370000000000000,0.380000000000000,0.390000000000000,0.400000000000000,0.410000000000000,0.420000000000000,0.430000000000000,0.440000000000000,0.450000000000000,0.460000000000000,0.470000000000000,0.480000000000000,0.490000000000000,0.500000000000000,0.510000000000000,0.520000000000000,0.530000000000000,0.540000000000000,0.550000000000000,0.560000000000000,0.570000000000000,0.580000000000000,0.590000000000000,0.600000000000000,0.610000000000000,0.620000000000000,0.630000000000000,0.640000000000000,0.650000000000000,0.660000000000000,0.670000000000000,0.680000000000000,0.690000000000000,0.700000000000000,0.710000000000000,0.720000000000000,0.730000000000000,0.740000000000000,0.750000000000000,0.760000000000000,0.770000000000000,0.780000000000000,0.790000000000000,0.800000000000000,0.810000000000000,0.820000000000000,0.830000000000000,0.840000000000000,0.850000000000000,0.860000000000000,0.870000000000000,0.880000000000000,0.890000000000000,0.900000000000000];
P_exp=[0.559256540754102,0.617268493908312,0.681907816813149,0.734998671894103,0.778495118965288,0.840514050227151,0.890770106313650,0.982184670626720,1.09253769506987,1.20290997765837,1.29554550838617,1.40324861203016,1.48489544502537,1.54546231212839,1.65291891151199,1.78419783668460,1.90919090321733,2.07926728799749,2.28597266521960,2.46774259585189,2.61924642526347,2.74742478903590,2.84764417739968,2.87158975531889,2.86745695732833,2.84935815233519,2.82144154484631,2.75232020958203,2.69405276503294,2.63990656358720,2.63548104178743,2.63082057061448,2.61753630195704,2.61945056160411,2.56473817068529,2.39869829316898,2.21621966278808,2.04257666949044,1.83503934038952,1.68035791699637,1.57032072613684,1.45544203753814,1.32972871636959,1.21777341423392,1.08868606600132,0.976904087173732,0.881672558453724,0.800938561547755,0.739324051338090,0.686225492998998,0.630064889550407,0.579870459529005,0.536535780878700,0.489195407997097,0.448033048141456,0.427723408062969,0.414373661711359,0.407818189036679,0.394094834665420,0.368323584567722,0.320212885872409,0.275468510983012,0.228008737600285,0.196972310565873,0.187115991779442,0.206481982736134,0.218664685479972,0.237113988718347,0.263620899968139,0.294310680386380,0.296240346549726,0.307028759570747,0.322200326471783,0.316711755049093,0.298924932010509,0.300357738024011,0.292993423244935,0.265196216257178,0.254091969652535,0.246296272417781,0.233285469424204,0.226067516549733,0.236686457891738,0.233343243860232,0.227288482964465,0.218537581720709,0.201967873467788,0.178980067310299,0.135456658835634,0.0966801692085699,0];
dCWBias=0.005;CWBias=0.001:dCWBias:CWBias_exp(end);
P1=interp1(CWBias_exp,P_exp,CWBias,'spline');
P_exp=[1.23407523693263,2.28235067261608,3.16814072439635,4.02097757124949,4.78197831207942,5.48006838895755,6.09821498249095,6.53950534523820,6.86466461718851,6.95005958849459,6.76100463197052,6.33813881989921,5.75439769142952,5.01836899383535,4.28215352149304,3.58483778159158,2.93880727461703,2.38874785906202,1.92177207730486,1.48926791000138,1.15077759040255,0.883483470110528,0.652781635891708,0.483871114182710,0.398880821497601,0.324602838846467,0.274072487193284,0.252258753066705,0.243056205580198,0.205549501469230,0.184284418165314,0.168058361920606,0.152081338673419,0.147186283815903,0.164050487116758,0.171035084469096,0.182183202248733,0.201849026771696,0.193234041513714,0.175358141660432,0.165918234598163,0.154820701646147,0.136034274895680,0.139404002643381,0.148015096760777,0.158649583983028,0.171809421445762,0.188759229839514,0.186113254240857,0.182031447765869,0.173074042136297,0.154454934431038,0.135730765929951,0.133633441153957,0.134559532613487,0.131901883593071,0.138279463013952,0.139516845720383,0.132773499084393,0.129170302901515,0.130571113512569,0.121851067458759,0.115430685491428,0.108344918483847,0.0879242126871503,0.0675035068904534,0.0539584465096793,0.0419270398169607,0.0353198831014898,0.0422305487826891,0.0527521929279384,0.0699938368656603,0.0844455330030332,0.0975080919511109,0.102037379593518,0.0989945076550625,0.0869825566652752,0.0785154347495716,0.0754141957023217,0.0761184921484349,0.0799707213288332,0.0856245486006702,0.0947375998536931,0.0929048726375643,0.0908231124239147,0.0894728866404822,0.0856012017571526,0.0740639699188892,0.0628341381869403,0.0489440633408512,0.0263624774719173];
% P = lognpdf(CWBias,log(0.318),0.298)/sum(lognpdf(CWBias,log(0.318),0.298))/dCWBias;
% P1 = lognpdf(CWBias,log(0.17),0.65)/sum(lognpdf(CWBias,log(0.17),0.65))/dCWBias;
P2 = interp1(CWBias_exp,P_exp,CWBias,'spline');
%  PWT= [lognpdf(CWBias,log(0.318),0.298)/sum(lognpdf(CWBias,log(0.318),0.298))/dCWBias lognpdf(CWBias,log(0.318),0.298)/sum(lognpdf(CWBias,log(0.318),0.298))/dCWBias]/2;
CWBias=[CWBias CWBias];
P=[P1 P2]/2;
%      CWBias=[0.2 0.2];P=[0.5 0.5];dCWBias=1;
% CC=winter(length(CWBias));
dt=0.1;    %s
dx=20;      %um

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
u=1;
A=500;
kA=280;
G0=9.3*0.84/60;% 10uM/min for OD 1 cells.
GSfactor=0.18;
GS=G0*GSfactor;
betafactor=1;
beta0=9.3*0.84/60/6;
beta=beta0*betafactor;
S=0;
dt=0.5;
T=round(120*60/dt);
for t=1:T
    u1=u;A1=A;A1(A1<=0)=0;
    S1=S; S1(S1<=0)=0;
    Gamma1=G0*u1./(1+(0.5./A1));%.*((kaO./(1+((Oex-S1)/ktO).^5))+1-kaO);
    GammaS=GS*u1./(1+(0.5./A1)); Gamma_beta=beta.*u1./(1+(0.5./S1))./(1+A1./kA);%Gamma_beta=beta.*u1;
      
    Aa=-Gamma1;
    Sa=GammaS-Gamma_beta;
    
    %%
    A1=A+Aa*dt/2;A1(A1<=0)=0;
   
    Gamma1=G0*u1./(1+(0.5./A1));%.*((kaO./(1+((Oex-S1)/ktO).^5))+1-kaO);
    S1=S+Sa*dt/2;S1(S1<=0)=0;
    GammaS=GS*u1./(1+(0.5./A1)); Gamma_beta=beta.*u1./(1+(0.5./S1))./(1+A1./kA);%Gamma_beta=beta.*u1;


    %f=ki'*A1;
    
    Ab=-Gamma1;
    Sb=GammaS-Gamma_beta;
   
    %%
   A1=A+Ab*dt/2;A1(A1<=0)=0;
    %Gamma1=G0*u1./(1+(0.5./A1)).*((kaO./(1+((Oex-S1)/ktO).^5))+1-kaO);
    Gamma1=G0*u1./(1+(0.5./A1));%.*((kaO./(1+((Oex-S1)/ktO).^5))+1-kaO);
    S1=S+Sb*dt/2;S1(S1<=0)=0;
    GammaS=GS*u1./(1+(0.5./A1)); Gamma_beta=beta.*u1./(1+(0.5./S1))./(1+A1./kA);%Gamma_beta=beta.*u1;
        %f=ki'*A1;
        Ac=-Gamma1;
    Sc=GammaS-Gamma_beta;
    
    
    %%
   A1=A+Ac*dt;A1(A1<=0)=0;
    %Gamma1=G0*u1./(1+(0.5./A1)).*((kaO./(1+((Oex-S1)/ktO).^5))+1-kaO);
    Gamma1=G0*u1./(1+(0.5./A1));%.*((kaO./(1+((Oex-S1)/ktO).^5))+1-kaO);
     S1=S+Sc*dt;S1(S1<=0)=0;
     GammaS=GS*u1./(1+(0.5./A1)); Gamma_beta=beta.*u1./(1+(0.5./S1))./(1+A1./kA);%Gamma_beta=beta.*u1;
       %f=ki'*A1;
     
    Ad=-Gamma1;
    Sd=GammaS-Gamma_beta;
   
     A=A+1/6*(Aa+2*Ab+2*Ac+Ad)*dt;
     S=S+1/6*(Sa+2*Sb+2*Sc+Sd)*dt;
    Af(:,t)=A;
  AS(:,t)=S;
  
    
end

plot((1:T)*dt,AS);
hold on;
plot((1-600:T-600)*dt,Af)
%%relantion TB vs mu and chi
dCWBias=0.005; CWBias=0.01:dCWBias:0.99;

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
Diff=v^2*(1-CWBias)./3./(lambdaR*theta+Drot);
lambdaRdiff=-(lambdaR).*a.*(a-1).*(alpha*g2*Kd)./(2*(Kd+alpha*a).^2); 
tau=20*exp(-6*CWBias);
N=6;
%  tau=0;
ki=theta*v^2*N*tau.*(lambdaRdiff./(1+lambdaR./lambdaT)./(lambdaR*theta+Drot))./(1+tau.*(lambdaR+Drot))/3;
plot(CWBias,Diff,'r',CWBias,ki,'g');
xlabel('Tumble bias'); ylabel('diff. coef. $$\mu$$, chemo. coef. $$\chi$$ ($$\mu^2$$/s)');
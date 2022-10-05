kaO=0.8;
ktO=80;
Oex=250;
O=0:250;

gamma=(kaO./(1+((Oex-O)/ktO).^5))+1-kaO;
plot(O,gamma);hold on;
ktO=100;
kappa=4;kaO=0.1;
gamma2=(kaO+(O/ktO).^kappa)./(1+(O/ktO).^kappa);
plot(O,gamma2);
% This code is for checking out single trajectories given a set of TB for two phenotypes.
% We can plot bacteria density (u) , Substrate (A) and gradient vs space and time, postion and intensity of
% peak of the wave vs time, the number density in the peak vs time with a linear consumption
% function.
%%
clear all;
close all;
%  dCWBias=0.1; CWBias=0.01:dCWBias:0.99;P = lognpdf(CWBias,log(0.3),0.3)/sum(lognpdf(CWBias,log(0.3),0.3))/dCWBias;

CWBias=[0.1 0.1]; % TB for two phenotypes
P=[0.5 0.5];      % Prop for two phenotypes
% Grid size
dt=0.05;    %s
dx=10;      %um
%plot(CWBias,P)

%% calculate the chemotactic coefficients based on the PLOS Comp.
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
Diff=v^2*(1-CWBias).*tauR/3*(1-0.1)/(1);
lambdaRdiff=-(lambdaR).*a.*(a-1).*(alpha*g2*Kd)./(2*(Kd+alpha*a).^2);
tau=20.8*exp(-5.86*CWBias)/2;
%tau=5;
N=6;
ki=v^2*N*tau.*(lambdaRdiff./(1+lambdaR./lambdaT)./(lambdaR*(1-0.1)+Drot))./(1+tau.*(lambdaR+Drot))/3*(1-0.1);
%%
Ki=3;      %uM
Ka=1000;    %uM
Asp=200;   %initial Asp conc in units of uM
L=20000/dx; %position
LX=0:dx:(L*dx-dx);
u=zeros(round(L),length(CWBias)); % bacteria variable
Drho=ones(round(L),1)*Diff/dx^2;  % diffusion constant for bacteria
A=ones(round(L),1)*Asp; %Substrate variable
T=200*60/dt; % time

% diffusion constant for substrate
%Dn=500/dx^2;
Dn=0;

% Consumption function
%G0=1/6.8;% 10uM/min for OD 1 cells.
%G0= log(2)/10;
G0=(log(2)/500);

%Parameters for initial distribution of bacteria
DL=1500/dx; V=600/dx; % Width for initial distribution for bacteria
ux=zeros(); % inital bacteria value
tot=zeros(); % variable to store number of cells in a peak
for x=1:round(DL)
    ux(x)=exp(-(x)^2/V^2); % initial distribution of bacteria
end
u(2:round(DL)+1,:)=ux'*P*500;   % populate initial bacteria value
u(1,:)=u(2,:);
%u(2:round(DL)+1,:)=ux'*P*500;%*dCWBias; %unit is OD1,~5*dx cells in
%600um*10um*dx
sum(sum(u(2:end-1,:))) % total density to check the conservation of cells. Note: for number of cells, multiply by dx

% Parameters for the peak
index=nan(1,2); % index position values for the peak
value= nan(1,2); % max value of  the peak

ua=u;ub=u;uc=u;ud=u;Aa=A;Ab=A;Ac=A;Ad=A;
%%
onemaxtrix=ones(1,length(CWBias));
% Empty matrices to plot
h0=[];h1=[];h2=[];h3=[];h4=[];h5=[];
tempT=zeros();
pos=nan(T,1); % Positon of total
pos1=nan(T,1); % position of first pheno
pos2=nan(T,1); %position of second pheno
TotU=nan(T,1); % total density of total
TotU1=nan(T,1); % total density of first pheno
TotU2=nan(T,1); % total density of second pheno
i=0;
g=u;
A1=A;

% for making a gif
%figure(1), hold on;
figure
filename = 'TB1_simple1.gif';


for t=1:T
    % Doing R-K integration scheme
    
    u1=u;A1=A;
    
    %Nonlinear consumption function
    %Gamma1=G0*sum(u1')'.*(0+1./(1+(sum(u1')'/10).^3))./(1+1./A1);
    %Gamma1=G0*sum(u1')'./(1+0.5./A1);
    %Gamma1=G0*sum(u1,2).*(1+(0.2/20)*sum(u1,2))./(1+(1/20).*sum(u1,2))./(1+0.5./A1);
    
    A1(A1<0)=0;                      %  Disallowing non zero values of A
    Gamma1=G0.*(A1>0.001).*sum(u1')'; % linear consumption function with threshold
    
    
    f=log(A1)*ki;          %linear receptor
%     f(f<0)=0;              %  Disallowing non zero values of f
    %f=log((1+A1/Ki)./(1+A1/Ka))*ki;   % nonlinear receptor
    
    
    g(2:end-1,:)=u1(2:end-1,:).*(f(3:end,:)-f(1:end-2,:))/2/dx;   g(1,:)=-g(2,:);g(end,:)=-g(end-1,:);
    ua(2:end-1,:)=-u1(2:end-1,:).*(2*Drho(2:end-1,:))+Drho(2:end-1,:).*(u1(1:end-2,:)+u1(3:end,:))-(g(3:end,:)-g(1:end-2,:))/2/dx;
    Aa(2:end-1)=-A1(2:end-1)*(2*Dn)+Dn*(A1(1:end-2)+A1(3:end))-Gamma1(2:end-1);
    
    Aa(1)=Aa(2); Aa(end)=Aa(end-1);
    ua(1,:)=ua(2,:);ua(end,:)=ua(end-1,:);
    
    %%
    u1=u+ua*dt/2;A1=A+Aa*dt/2;
    A1(A1<0)=0;
    Gamma1=G0.*(A1>0.001).*sum(u1')';
    
    
    f=log(A1)*ki;
%     f(f<0)=0;
    %f=log((1+A1/Ki)./(1+A1/Ka))*ki;
    g(2:end-1,:)=u1(2:end-1,:).*(f(3:end,:)-f(1:end-2,:))/2/dx;g(1,:)=-g(2,:);g(end,:)=-g(end-1,:);
    ub(2:end-1,:)=-u1(2:end-1,:).*(2*Drho(2:end-1,:))+Drho(2:end-1,:).*(u1(1:end-2,:)+u1(3:end,:))-(g(3:end,:)-g(1:end-2,:))/2/dx;
    Ab(2:end-1)=-A1(2:end-1)*(2*Dn)+Dn*(A1(1:end-2)+A1(3:end))-Gamma1(2:end-1);
    
    Ab(1)=Ab(2); Ab(end)=Ab(end-1);
    ub(1,:)=ub(2,:);ub(end,:)=ub(end-1,:);
    %%
    u1=u+ub*dt/2;A1=A+Ab*dt/2;
    A1(A1<0)=0;
    Gamma1=G0.*(A1>0.001).*sum(u1')';
    
    
    f=log(A1)*ki;
%     f(f<0)=0;
    %    f=log((1+A1/Ki)./(1+A1/Ka))*ki;
    g(2:end-1,:)=u1(2:end-1,:).*(f(3:end,:)-f(1:end-2,:))/2/dx;g(1,:)=-g(2,:);g(end,:)=-g(end-1,:);
    uc(2:end-1,:)=-u1(2:end-1,:).*(2*Drho(2:end-1,:))+Drho(2:end-1,:).*(u1(1:end-2,:)+u1(3:end,:))-(g(3:end,:)-g(1:end-2,:))/2/dx;
    Ac(2:end-1)=-A1(2:end-1)*(2*Dn)+Dn*(A1(1:end-2)+A1(3:end))-Gamma1(2:end-1);
    Ac(1)=Ac(2);
    Ac(end)=Ac(end-1);
    uc(1,:)=uc(2,:);uc(end,:)=uc(end-1,:);
    
    %%
    u1=u+uc*dt;A1=A+Ac*dt;
    A1(A1<0)=0;
    Gamma1=G0.*(A1>0.001).*sum(u1')';
    
    
    f=log(A1)*ki;
    f(f<0)=0;
    %    f=log((1+A1/Ki)./(1+A1/Ka))*ki;
    
    g(2:end-1,:)=u1(2:end-1,:).*(f(3:end,:)-f(1:end-2,:))/2/dx;g(1,:)=-g(2,:);g(end,:)=-g(end-1,:);
    ud(2:end-1,:)=-u1(2:end-1,:).*(2*Drho(2:end-1,:))+Drho(2:end-1,:).*(u1(1:end-2,:)+u1(3:end,:))-(g(3:end,:)-g(1:end-2,:))/2/dx;
    Ad(2:end-1)=-A1(2:end-1)*(2*Dn)+Dn*(A1(1:end-2)+A1(3:end))-Gamma1(2:end-1);
    Ad(1)=Ad(2); Ad(end)=Ad(end-1);
    ud(1,:)=ud(2,:);ud(end,:)=ud(end-1,:);
    %%
    
    u=u+1/6*(ua+2*ub+2*uc+ud)*dt;
    A=A+1/6*(Aa+2*Ab+2*Ac+Ad)*dt;
    A(1)=A(2); A(end)=A(end-1);g(1,:)= -g(2,:);g(end,:)=-g(end-1,:);
    u(1,:)=u(2,:);u(end,:)=u(end-1,:);
    A(A<0)=0;
    
    %Plotting every 10 mins
    
    if (mod(t,600/dt)==0 && t> 6000)
        i=i+1;
        %         if ~isempty(h1) && ~isempty(h2) && ~isempty(h3)&& ~isempty(h4)
        %             delete(h1);delete(h2);delete(h3);delete (h4);
        %         end
        
        if  ~isempty(h3)&& ~isempty(h4)
            delete(h3);delete(h4);
        end
        
        %plot u, A profiles
        subplot(2,3,1)
        %h1=plot(LX,u(:,1),'.-','MarkerSize',1,'Color',[0 0.7 0]); hold on;
        %h2= plot(LX,u(:,2),'.-','MarkerSize',1,'Color',[0 0 1]); hold on;
        %h2= plot(diff(LX),diff(A),'.-','MarkerSize',1,'Color',[0 0 1]); hold on;
        h4 = plot(LX,sum(u,2),'.--','MarkerSize',5,'Color',[0.5 0 0]); hold on;
        h3=plot(LX,A,'r');
        
        
        
        %finding maxima for the peak to plot peak position later
        [pks1, loc1]=findpeaks(u(:,1));
        [pks2, loc2]=findpeaks(u(:,2));
        [pks, loc]= findpeaks(sum(u,2));
        
        
        % If multiple peak posotion are returned, then pick the last one as
        % that is the actual wave; if no peak posotions are returned then
        % origin in the 'peak' position
        if ~isempty(loc1)
            index(:,1)=loc1(end);
            value(:,1)= pks1(end);
        else
            index(:,1)=1;
            value(:,1)= max(u(:,1));
        end
        
        if ~isempty(loc2)
            index(:,2)=loc2(end);
            value(:,2)= pks2(end);
        else
            index(:,2)=1;
            value(:,2)= max(u(:,2));
        end
        
        if ~isempty(loc)
            ind3=loc(end);
            val3= pks(end);
        else
            ind3=1;
            val3= max(sum(u,2));
        end
        
        xlim([0 L*dx]);ylim([0 500]);
        xlabel('x(microns)');
        ylabel('u, A');
        title('Time evolution of u and A');
        %legend({'LowTB','HighTB','Total'},'Location','northeast','FontSize',6);
        legend({'Bacteria','Substrate'},'Location','northeast','FontSize',6);
        axis square;
        box on;
        drawnow;
        
        
        
        
        % Storing the peak posistion values
        tempT(i)=t;
        pos1(t)= index(:,1)*dx;
        pos2(t)= index(:,2)*dx;
        pos(t)= ind3*dx;
        
        %Plot peak position vs time by using maxima for the peak
        subplot(2,3,2)
        %   plot(t*dt,index(:,1)*dx,'.','MarkerSize',5,'Color',[0,0.7,0]),hold on
        %    plot(t*dt,index(:,2)*dx,'.','MarkerSize',5,'Color',[0,0,1]),hold on
        plot(t*dt,ind3*dx,'.','MarkerSize',5,'Color',[0.5,0,0]),hold on;
        xlabel('t(sec)');
        ylabel('x(microns)');
        title('Peak position vs Time ')
        %  legend({'LowTB','HighTB','Total'},'Location','northeast','FontSize',6);
        axis square;
        box on
        drawnow;
        
        %Plot peak intensity vs time by finding maxima for the peak
        subplot(2,3,3)
        %    h2= plot(diff(A),'.-','MarkerSize',1,'Color',[0 0 1]); hold on;
        %    plot(t*dt,value(:,1),'.','MarkerSize',5,'Color',[0,0.7,0]),hold on;
        %  plot(t*dt,value(:,2),'.','MarkerSize',5,'Color',[0,0,1]),hold on;
        plot(t*dt,val3,'.','MarkerSize',5,'Color',[0.5,0,0]),hold on;
        xlabel('t');
        ylabel('Max u');
        title('Peak intensity vs Time ')
        %  legend({'LowTB','HighTB','Total'},'Location','northeast','FontSize',6);
        axis square;
        box on
        drawnow;
        
        if ~isempty(h5)
            delete(h5);
        end
        
        %Plot gradient vs time
        subplot(2,3,4)
        h5= plot(diff(A)./(dx*(A(1:end-1)+A(2:end)/2)),'.-','MarkerSize',1,'Color',[0 0 1]); hold on;
        xlabel('x');
        ylabel('diff(A)/A');
        title('diff(A)/A vs x ')
        xlim([0 L]);ylim([0 0.01]);
        % legend({'LowTB','HighTB','Total'},'Location','northeast','FontSize',6);
        axis square;
        box on
        drawnow;
        
        
        
        
        % The next section is to figure out how many cells are in the wave.
        % One has to determine how to 'choose' what is the wave. If there
        % is a single peak, then density in the wave is the total number of 
        % cells. In case of multiple waves, finding the number of cells in 
        % between the last two minimas would be the way to go. 
        
        Data=sum(u,2);
        DataInv = 1.01*max(Data) - Data;
        [Min,MinIdx] = findpeaks(DataInv);
        
        % for a single wave; total number of cells
        TotU(t)=sum(Data(2:end-1));
        TotU1(t)=sum(u(2:end-1,1));
        TotU2(t)=sum(u(2:end-1,2));
        
        % for multiple wave; the number of cells within the last two minima
        %         if ~isempty(MinIdx)
        %             Min = Data(MinIdx(1));
        %             TotU(t)=sum(Data(MinIdx(1):end-1));
        %             TotU1(t)=sum(u(MinIdx(1):end-1,1));
        %             TotU2(t)=sum(u(MinIdx(1):end-1,2));
        %         else
        %             TotU(t)=0;
        %             TotU1(t)=0;
        %             TotU2(t)=0;
        %             TotU(t)=sum(Data(1:end));
        %             TotU1(t)=sum(u(1:end,1));
        %             TotU2(t)=sum(u(1:end,2));
        %             TotU(t)=sum(Data(2:end-1));
        %             TotU1(t)=sum(u(2:end-1,1));
        %             TotU2(t)=sum(u(2:end-1,2));
        %         end
        
        
        %Plotting the total density for cells in the wave
        subplot(2,3,5)
        %plot(t*dt,TotU1(t),'.','MarkerSize',5,'Color',[0,0.7,0]),hold on;
        %plot(t*dt,TotU2(t),'.','MarkerSize',5,'Color',[0,0,1]),hold on;
        plot(t*dt,TotU(t),'.','MarkerSize',5,'Color',[0.5,0,0]),hold on;
        xlabel('t');
        ylabel('Bacteria density');
        title('Peak Density vs Time ')
        % legend({'LowTB','HighTB','Total'},'Location','northeast','FontSize',6);
        axis square;
        box on
        drawnow;
        %
        
        subplot(2,3,6)
        % to plot peak intensity vs time
        % plot(t*dt,value(:,1),'.','MarkerSize',5,'Color',[0,0.7,0]),hold on;
        % plot(t*dt,value(:,2),'.','MarkerSize',5,'Color',[0,0,1]),hold on;
        % semilogy(t*dt,val3,'.','MarkerSize',5,'Color',[0.5,0,0]),hold on;
        % xlabel('t');
        % ylabel('Log (Max u)');
        % title('Peak intensity vs Time ')
        % legend({'LowTB','HighTB','Total'},'Location','northeast','FontSize',6);
        % axis square;
        % box on
        % drawnow;
        
        % Plotting bacteria density, gradient, substrate on the same axis.
        plot(u/max(u))
        hold on; plot(diff(A)./A(2:end)/max(diff(A)./A(2:end)))
        plot (A/max(A))
        xlabel('t');
        ylabel('Norm A, g, u');
        title('Norm A, g, u ')
        %  legend({'LowTB','HighTB','Total'},'Location','northeast','FontSize',6);
        axis square;
        box on
        drawnow;hold off;
        
        % Make the gifs
        frame = getframe(gcf);
        im = frame2im(frame);
        [Amp,map] = rgb2ind(im,256);
        if i == 1;
            imwrite(Amp,map,filename,'gif','LoopCount',Inf,'DelayTime',0.2);
        else
            imwrite(Amp,map,filename,'gif','WriteMode','append','DelayTime',0.2);
        end
    end
end

ind=~isnan(pos);
pos=pos(ind);

ind1=~isnan(pos1);
pos1=pos1(ind1);

ind2=~isnan(pos2);
pos2=pos2(ind2);

tempT= tempT(:);

% Plot peak position vs time to calculate velocity later
% figure(2)
% subplot(1,2,1)
% plot(tempT*dt,pos1,'.-g'), hold on;
% plot(tempT*dt,pos2,'.-b');
% %plot(tempT*dt,pos,'.-r');
% xlabel('t(sec)')
% ylabel('Position(microns)');
% title('Peak position vs Time ')
% legend({'LowTB','HighTB','Total'},'Location','northwest','FontSize',6);
% axis square;box on;


% calculating velocity
m1=diff(pos1); % change in Bacteria 1
m3= diff(pos2); % change in Bacteria 2
m4= diff(pos); % change in Total Bacteria
m2= diff(tempT)*dt; % change in time
p = (tempT(1:end-1)+tempT(2:end))*dt/2; % Velocity= dx/dt

% Plot velcoty of the wave vs time
% subplot(1,2,2)
% plot(p,m1./m2,'.-g'), hold on;
% plot(p,m3./m2,'.-b'),hold on;
% %plot(p,m4./m2,'.-r'),
% xlim([min(p) max(p)]);ylim([0 5]);
% xlabel('t(sec)')
% ylabel('Velocity(microns/s)');
% title('Peak velocity vs Time ')
% legend({'LowTB','HighTB','Total'},'Location','northeast','FontSize',6);
% axis square;box on;

% To check if the density is conserved
sum(sum(u(2:end-1,:)))


ind=~isnan(TotU);
TotU=TotU(ind);

ind1=~isnan(TotU1);
TotU1=TotU1(ind1);

ind2=~isnan(TotU2);
TotU2=TotU2(ind2);

% Plot bacteria density vs time
% figure(3)
% plot(tempT*dt,TotU1,'.-g'), hold on;
% plot(tempT*dt,TotU2,'.-b');
% plot(tempT*dt,TotU,'.-r');
% xlabel('t(sec)')
% ylabel('Density(microns)');
% title('Peak Density vs Time ')
% legend({'LowTB','HighTB','Total'},'Location','northwest','FontSize',6);
% axis square;box on;


% Wave speed from classic K-S for a constant consumption function
mean(TotU)*dx*G0/(Asp);





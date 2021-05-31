clear all
close all

%---------------%
%-D parameters--%
%---------------%
param.L=500;    %--normal fault length (m)
param.gamma=60; %--normal fault angle (deg)
param.z=3e3;    %--depth (m)
param.rho=3e3;  %--rock density (kg/m3)
param.g=9.81;   %--gravity (m2/s)
param.mu=3e10;  %--shear modulus (Pa)
param.f0=0.6;   %--static friction coeff.
param.a=0.001;  %--friction param. 1 (rate-and-state)
param.b=0.002;  %--friction param. 2 (rate-and-state)
param.dc=0.1e-3;%--critical slip (m)
param.vp=1e-9;  %--tectonic extension rate (m/s)
param.alpha=0.3;%--Linker-Dieterich coeff.
param.u0=-1.0;  %--initial displacement (controls initial (compressive) stress)
param.wd=0.1;   %--dyke width (m)
param.dtau1=0.0*param.rho*param.g*param.z; %--ça je me souviens plus ce que c'est
param.taud=0.5*param.rho*param.g*param.z;  %--extensional fracture criterion (for dyke amplacement)
param.niter=100000; %--number of iterations

format long

year=3600*24*365;  %--number of seconds in one year

%---------------------%
%-Non Dim parameters--%
%---------------------%
par.r=param.a/param.b;
par.b=param.b;
par.f0=param.f0/param.b;
%par.beta=0;
par.beta=0.5*sqrt(param.mu*param.rho)*param.vp/(param.b*param.rho*param.g*param.z); %--normalized radiation damping
par.k=param.mu*param.dc/(param.L*param.b*param.rho*param.g*param.z);
par.ks=param.mu*param.dc/(param.L*param.rho*param.g*param.z);
par.wd=param.wd/param.dc;
par.taud=param.taud/(param.rho*param.g*param.z);
par.dtau1=param.dtau1/(param.rho*param.g*param.z);
par.c=cos(param.gamma*pi/180);
par.s=sin(param.gamma*pi/180);
par.alpha=param.alpha;
par.u0=param.u0/param.dc;


disp(['Normalized stiffness : ',num2str(par.k),' Normalized critical stiffness : ',num2str(1-par.r)]);

save parametres.mat param par

%--------------------%
%-Comput parameters--%
%--------------------%
par.err0=1e-10;
par.dtmin=1e-10;
par.dtmax=1e2;
par.maxiterrk=10;
par.maxiternew=30000;
par.safe=0.8;
repres='./results/';


%---------------------%
%-Initial conditions--%
%---------------------%
t=0;            %-Norm. time tvp/dc
dt=1e-6;        %-Norm. Time step
th=log(0.001);  %-Norm log state var
u=0;            %-Norm. slip u/dc
sig=1-par.k*par.u0*par.s;  %-Norm. normal stress sig/(rho*g*z)
w=0.0;          %--Cumulative dyke width
dw=0.0;         %--Rate of cumulative dyke width increase
us=0.0; 
tau1=0.0;
phi=-th;        %-Initial guess for phi (log normalized slip rate)
phi=phinew(par,phi,th,u,w,tau1,t);  %-compute initial slip rate (from stress balance)


fco=friction(par,phi,th);  %friction coefficient
tau=fco*sig;    %-Norm. shear stress tau/(rho*g*z)

%%

%---------------------%
%-Open output files---%
%---------------------%
system(['rm -rf ',repres]);
system(['mkdir ',repres]);
fidt=fopen([repres,'t.data'],'w+');
fidv=fopen([repres,'v.data'],'w+');
fidth=fopen([repres,'th.data'],'w+');
fidu=fopen([repres,'u.data'],'w+');
fidw=fopen([repres,'w.data'],'w+');
fids=fopen([repres,'s.data'],'w+');
fidtau=fopen([repres,'tau.data'],'w+');

for iter=1:param.niter
    
    %----------------------%
    %-Write output files---%
    %----------------------%
    fwrite(fidt,t,'real*8');
    fwrite(fidv,exp(phi),'real*8');
    fwrite(fidth,exp(th),'real*8');
    fwrite(fidu,u,'real*8');
    fwrite(fidw,w,'real*8');
    fwrite(fids,sig,'real*8');
    fwrite(fidtau,tau,'real*8');
    %----------------------%
    %-Update var.----------%
    %----------------------%
    phiex=phi;thex=th;uex=u;wex=w;dwex=dw;sigex=sig;tauex=tau;
    [phi,th,u,w,dw,sig,dt,error]=rk4fehlberg(dt,par,phi,th,u,w,dw,sig,tau1,t);
    phi=phinew(par,phi,th,u,w,tau1,t);
    fco=friction(par,phi,th);
    tau=fco*sig;
    
    %------------------------%
    %--Check dike criterion--%
    %------------------------%
    if ((par.ks*(u-us)>=par.taud) & (par.ks*(uex-us)<par.taud))
        
        par.td=dt*((par.taud/par.ks-(uex-us))/(u-uex))+t;
        
        phi=phiex+(phi-phiex)*(par.td-t)/dt;
        th=thex+(th-thex)*(par.td-t)/dt;
        u=uex+(u-uex)*(par.td-t)/dt;
        w=wex+(w-wex)*(par.td-t)/dt;
        dw=dwex+(dw-dwex)*(par.td-t)/dt;
        sig=sigex+(sig-sigex)*(par.td-t)/dt;
        tau=tauex+(tau-tauex)*(par.td-t)/dt;
        t=par.td;
        us=u;
        
        
        %----------------------%
        %-Write output files---%
        %----------------------%
        fwrite(fidt,t,'real*8');
        fwrite(fidv,exp(phi),'real*8');
        fwrite(fidth,exp(th),'real*8');
        fwrite(fidu,u,'real*8');
        fwrite(fidw,w,'real*8');
        fwrite(fids,sig,'real*8');
        fwrite(fidtau,tau,'real*8');
        
        disp('Oh le joli dyke');
        %par.td=par.td+par.T;
        w=w+par.wd;
        tau1=tau1+par.dtau1;
        th=th-(par.alpha/par.b)*log(1+par.b*par.k*par.s*par.wd/sig);
        sig=sig+par.b*par.k*par.s*par.wd;
        phi=phinew(par,phi,th,u,w,tau1,t);
        fco=friction(par,phi,th);
        tau=fco*sig;
    else
        t=t+dt;
    end
    
    if (mod(floor(100*iter/param.niter),10)==0 & mod(floor(100*(iter-1)/param.niter),10)>0)
        disp([num2str(floor(100*iter/param.niter)),' % completed']);
    end
    
    
end
%----------------------%
%-Close output files---%
%----------------------%
fclose(fidt);
fclose(fidv);
fclose(fidth);
fclose(fidu);
fclose(fidw);
fclose(fids);
fclose(fidtau);


function [phi,th,u,w,dw,sig,dt,error]=rk4fehlberg(dt,par,phi,th,u,w,dw,sig,tau1,t)

passtep=0;
iterrk=0;

while ((passtep==0) & (iterrk<=par.maxiterrk))
    
    iterrk=iterrk+1;

    l1=grns(par,u,phi,th,w,dw,sig,tau1,t);
    m1=hrns(par,u,phi,th,w,tau1,t);
    
    w1=dykeopen(t+0.25*dt,w,par);
    dw1=dykeopenrate(t+0.25*dt);
    sig1=sigma(par,u+0.25*dt*m1,w1,t+0.25*dt);
    
    l2=grns(par,u+0.25*dt*m1,phi,th+0.25*dt*l1,w1,dw1,sig1,tau1,t+0.25*dt);
    m2=hrns(par,u+0.25*dt*m1,phi,th+0.25*dt*l1,w1,tau1,t+0.25*dt);
    
    w2=dykeopen(t+(3/8)*dt,w,par);
    dw2=dykeopenrate(t+(3/8)*dt);
    sig2=sigma(par,u+(3/32)*dt*m1+(9/32)*dt*m2,w2,t+(3/8)*dt);
    
    l3=grns(par,u+(3/32)*dt*m1+(9/32)*dt*m2,phi,th+(3/32)*dt*l1+(9/32)*dt*l2,w2,dw2,sig2,tau1,t+(3/8)*dt);
    m3=hrns(par,u+(3/32)*dt*m1+(9/32)*dt*m2,phi,th+(3/32)*dt*l1+(9/32)*dt*l2,w2,tau1,t+(3/8)*dt);
    
    w3=dykeopen(t+(12/13)*dt,w,par);
    dw3=dykeopenrate(t+(12/13)*dt);
    sig3=sigma(par,u+(1932/2197)*dt*m1-(7200/2197)*dt*m2+(7296/2197)*dt*m3,w3,t+(12/13)*dt);
    
    l4=grns(par,u+(1932/2197)*dt*m1-(7200/2197)*dt*m2+(7296/2197)*dt*m3,phi,th+(1932/2197)*dt*l1-(7200/2197)*dt*l2+(7296/2197)*dt*l3,w3,dw3,sig3,tau1,t+(12/13)*dt);
    m4=hrns(par,u+(1932/2197)*dt*m1-(7200/2197)*dt*m2+(7296/2197)*dt*m3,phi,th+(1932/2197)*dt*l1-(7200/2197)*dt*l2+(7296/2197)*dt*l3,w3,tau1,t+(12/13)*dt);
    
    w4=dykeopen(t+dt,w,par);
    dw4=dykeopenrate(t+dt);
    sig4=sigma(par,u+(439/216)*dt*m1-8*dt*m2+(3680/513)*dt*m3-(845/4104)*dt*m4,w4,t+dt);
    
    l5=grns(par,u+(439/216)*dt*m1-8*dt*m2+(3680/513)*dt*m3-(845/4104)*dt*m4,phi,th+(439/216)*dt*l1-8*dt*l2+(3680/513)*dt*l3-(845/4104)*dt*l4,w4,dw4,sig4,tau1,t+dt);
    m5=hrns(par,u+(439/216)*dt*m1-8*dt*m2+(3680/513)*dt*m3-(845/4104)*dt*m4,phi,th+(439/216)*dt*l1-8*dt*l2+(3680/513)*dt*l3-(845/4104)*dt*l4,w4,tau1,t+dt);
    
    w5=dykeopen(t+0.5*dt,w,par);
    dw5=dykeopenrate(t+0.5*dt);
    sig5=sigma(par,u-(8/27)*dt*m1+2*dt*m2-(3544/2565)*dt*m3+(1859/4104)*dt*m4-(11/40)*dt*m5,w5,t+0.5*dt);
    
    l6=grns(par,u-(8/27)*dt*m1+2*dt*m2-(3544/2565)*dt*m3+(1859/4104)*dt*m4-(11/40)*dt*m5,phi,th-(8/27)*dt*l1+2*dt*l2-(3544/2565)*dt*l3+(1859/4104)*dt*l4-(11/40)*dt*l5,w5,dw5,sig5,tau1,t+0.5*dt);
    m6=hrns(par,u-(8/27)*dt*m1+2*dt*m2-(3544/2565)*dt*m3+(1859/4104)*dt*m4-(11/40)*dt*m5,phi,th-(8/27)*dt*l1+2*dt*l2-(3544/2565)*dt*l3+(1859/4104)*dt*l4-(11/40)*dt*l5,w5,tau1,t+0.5*dt);
    

     
    verror=dt*abs((16/135-25/216)*[l1;m1]+(6656/12825-1408/2565)*[l3;m3]+(28561/56430-2197/4104)*[l4;m4]-(9/50-1/5)*[l5;m5]+(2/55)*[l6;m6]);
    error=max(verror);
    %(16/135-25/216)*[k1;l1;m1]+(6656/12825-1408/2565)*[k3;l3;m3]+(28561/56430-2197/4104)*[k4;l4;m4]-(9/50-1/5)*[k5;l5;m5]+(2/55)*[k6;l6;m6]
    %(16/135-25/216)*[k1;l1;m1]+(6656/12825-1408/2565)*[k3;l3;m3]+(28561/56430-2197/4104)*[k4;l4;m4]-(9/50-1/5)*[k5;l5;m5]+(2/55)*[k6;l6;m6]
    
    if error<=par.err0
        passtep=1;
        dt=min([par.safe*dt*((par.err0/error)^(1/4)) par.dtmax]);
    else
        dt=max([par.safe*dt*((par.err0/error)^(1/4)) par.dtmin]);
    end
    
    %verror
    
    
end

if iterrk==par.maxiterrk
    disp([' RK algo did not find solution to sufficient accuracy, error = ',num2str(error)]);
end

th=th+dt*((16/135)*l1+(6656/12825)*l3+(28561/56430)*l4-(9/50)*l5+(2/55)*l6);
u=u+dt*((16/135)*m1+(6656/12825)*m3+(28561/56430)*m4-(9/50)*m5+(2/55)*m6);
w=w4;
dw=dw4;
sig=sigma(par,u,w,t+dt);

end


function g=grns(par,u,phi,th,w,dw,sig,tau1,t)
phi=phinew(par,phi,th,u,w,tau1,t);
sig=sigma(par,u,w,t); 
g=exp(-th)-exp(phi)-par.alpha*par.k*par.b*(par.c*exp(phi)-1+dw)*par.s/sig;
end

function h=hrns(par,u,phi,th,w,tau1,t)
phi=phinew(par,phi,th,u,w,tau1,t);
h=exp(phi);
end

function sig=sigma(par,u,w,t) % contrainte normale
sig=1+par.k*(u*par.c-t+w-par.u0)*par.s; % w = épaisseur cumulée des dikes, à enlever ; ajouter ici terme sig_n
end

function phi=phinew(par,phi,th,u,w,tau1,t) % ajouter terme de shear stress (attention conditions)
if par.s==0
    errphi=1;iternew=0;
    while ((errphi>par.err0) & (iternew<par.maxiternew))
        iternew=iternew+1;
        phi0=phi;
        vphi(iternew)=phi;
        f=par.r*phi+par.beta*exp(phi)+par.f0+th+par.k*(par.c*u+w-t-par.u0)*par.c-tau1/par.b; 
        fp=par.r+par.beta*exp(phi);
        phi=phi-f/fp;
        errphi=abs(phi-phi0);
    end
    if iternew==par.maxiternew
        disp([' NR algo did not find solution to sufficient accuracy, error = ',num2str(errphi),' ',num2str(iternew),' iterations']);
        save vphi.mat vphi
        stop
    end
else
    errphi=1;iternew=0;
    M=1+par.b*par.k*(par.c*u+w-t-par.u0)*par.s;
    if M<0
        fm=par.r*M*log(-par.r*M/par.beta)+par.beta^2/(M*par.r)+M*(par.f0+th)-par.s/par.b+(M-1)*par.c/(par.b*par.s)-tau1/par.b;
        disp(['error : M < 0 ',num2str(M),' fm : ',num2str(fm)]);
        pause
    else
        while ((errphi>par.err0) & (iternew<par.maxiternew))
            iternew=iternew+1;
            phi0=phi;
            f=par.r*M*phi+par.beta*exp(phi)+M*(par.f0+th)-par.s/par.b+(M-1)*par.c/(par.b*par.s)-tau1/par.b;
            fp=par.r*M+par.beta*exp(phi);
            phi=phi-f/fp;
            errphi=abs(phi-phi0);
        end
        if iternew==par.maxiternew
            disp([' NR algo did not find solution to sufficient accuracy, error: ',num2str(errphi),' ',num2str(iternew),' iterations']);
            pause
        end
    end
end
end

function fco=friction(par,phi,th)
fco=par.b*(par.f0+par.r*phi+th);
end

function w=dykeopen(t,w,par)
w=w;
end

function dw=dykeopenrate(t)
dw=0;
end



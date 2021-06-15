%clear all
%close all
year=3600*24*365;  %--number of seconds in one year
%---------------%
%-D parameters--%
%---------------%
param.L=Lf;    %--normal fault length (m)
param.gamma=thetaf; %--normal fault angle (deg)
param.z=abs(xsi(2));    %--depth (m)
param.rho=3e3;  %--rock density (kg/m3)
param.g=9.81;   %--gravity (m2/s)
param.mu=mu;  %--shear modulus (Pa)
param.f0=0.6;   %--static friction coeff.
param.a=0.002;  %--friction param. 1 (rate-and-state)
param.b=0.001;  %--friction param. 2 (rate-and-state)
param.dc=0.1e-3;%--critical slip (m)
param.vp=1e-9;  %--tectonic extension rate (m/s)
param.alpha=0.3;%--Linker-Dieterich coeff.
param.u0=-0.5;  %--initial displacement (controls initial (compressive) stress)
param.dtau1=0.0*param.rho*param.g*param.z; %--ça je me souviens plus ce que c'est
param.niter=100000; %--number of iterations
param.tfinal=60*year; %--temps de modélisation en année
param.tdike=20*year;
param.nbd=1*year;

format long

%--------------------%
%-Termes de forçage--%
%--------------------%

%sigm.n=sigm.n/(param.rho*param.g*param.z);
%sigm.t=sigm.t/(param.rho*param.g*param.z);

sigm.n=0;
sigm.t=0;
 
%---------------------%
%-Non Dim parameters--%
%---------------------%
par.r=param.a/param.b;
par.b=param.b;
par.f0=param.f0/param.b;
par.beta=0.5*sqrt(param.mu*param.rho)*param.vp/(param.b*param.rho*param.g*param.z); %--normalized radiation damping
par.k=param.mu*param.dc/(param.L*param.b*param.rho*param.g*param.z);
par.ks=param.mu*param.dc/(param.L*param.rho*param.g*param.z);
par.dtau1=param.dtau1/(param.rho*param.g*param.z);
par.c=cos(param.gamma*pi/180);
par.s=sin(param.gamma*pi/180);
par.alpha=param.alpha;
par.u0=param.u0/param.dc;
par.tfinal=param.tfinal*param.vp/param.dc;
par.tdike=param.tdike*param.vp/param.dc;
par.nbd=param.nbd*param.vp/param.dc;
par.delayd=1*year*param.vp/param.dc;
par.wd=b/param.dc;
par.dsigmn=dsigm.n/(param.rho*param.g*param.z);
par.dsigmt=dsigm.t/(param.rho*param.g*param.z);


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
sig=1-par.k*par.u0*par.s+sigm.n;  %-Norm. normal stress sig/(rho*g*z)
us=0.0; 
tau1=0.0;
phi=-th;        %-Initial guess for phi (log normalized slip rate)
phi=phinew(par,phi,th,u,tau1,t,sigm);  %-compute initial slip rate (from stress balance)


fco=friction(par,phi,th);  %friction coefficient
tau=fco*sig;    %-Norm. shear stress tau/(rho*g*z)
tpsi=1;
tpsii=[1:10];

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
fids=fopen([repres,'s.data'],'w+');
fidtau=fopen([repres,'tau.data'],'w+');

while t<par.tfinal
    
    %----------------------%
    %-Write output files---%
    %----------------------%
    fwrite(fidt,t,'real*8');
    fwrite(fidv,exp(phi),'real*8');
    fwrite(fidth,exp(th),'real*8');
    fwrite(fidu,u,'real*8');
    fwrite(fids,sig,'real*8');
    fwrite(fidtau,tau,'real*8');
    %----------------------%
    %-Update var.----------%
    %----------------------%
    phiex=phi;thex=th;uex=u;sigex=sig;tauex=tau;
    [phi,th,u,sig,dt,error]=rk4fehlberg(dt,par,phi,th,u,sig,tau1,t);
    [dbi,dbc]=evoldyke(par,t);
    sigm.n=dbc*par.dsigmn;
    sigm.t=dbc*par.dsigmt;
    phi=phinew(par,phi,th,u,tau1,t,sigm);
    %disp(sig)
    fco=friction(par,phi,th);
    tau=fco*sig;

    t=t+dt;
    %sigm.n=sigm.n+dsigm.n/(param.rho*param.g*param.z)*dt*param.dc/param.vp;
    %sigm.t=sigm.t+dsigm.t/(param.rho*param.g*param.z)*dt*param.dc/param.vp;
    %sigm.n=sigm.n+dsigm.n/(param.rho*param.g*param.z)*b*dbi*dt*param.dc/param.vp;
    %sigm.t=sigm.t+dsigm.t/(param.rho*param.g*param.z)*b*dbi*dt*param.dc/param.vp;
    %if (mod(floor(100*iter/param.niter),10)==0 & mod(floor(100*(iter-1)/param.niter),10)>0)
    %    disp([num2str(floor(100*iter/param.niter)),' % completed']);
    %end
    %[dbi,dbc]=evoldyke(par,t);
    %sigm.n=dbc*par.dsigmn;
    %sigm.t=dbc*par.dsigmt;
    %disp(dbc)
    
    if t>tpsii(tpsi)/10*par.tfinal
            disp([num2str(floor(100*t/par.tfinal)),' % completed']);
            tpsi=tpsi+1;
    end
    
    
end
%----------------------%
%-Close output files---%
%----------------------%
fclose(fidt);
fclose(fidv);
fclose(fidth);
fclose(fidu);
fclose(fids);
fclose(fidtau);


function [phi,th,u,sig,dt,error]=rk4fehlberg(dt,par,phi,th,u,sig,tau1,t)

passtep=0;
iterrk=0;

while ((passtep==0) & (iterrk<=par.maxiterrk))
    
    iterrk=iterrk+1;

    l1=grns(par,u,phi,th,sig,tau1,t);
    m1=hrns(par,u,phi,th,tau1,t);
    
    sig1=sigma(par,u+0.25*dt*m1,t+0.25*dt);
    
    l2=grns(par,u+0.25*dt*m1,phi,th+0.25*dt*l1,sig1,tau1,t+0.25*dt);
    m2=hrns(par,u+0.25*dt*m1,phi,th+0.25*dt*l1,tau1,t+0.25*dt);

    sig2=sigma(par,u+(3/32)*dt*m1+(9/32)*dt*m2,t+(3/8)*dt);
    
    l3=grns(par,u+(3/32)*dt*m1+(9/32)*dt*m2,phi,th+(3/32)*dt*l1+(9/32)*dt*l2,sig2,tau1,t+(3/8)*dt);
    m3=hrns(par,u+(3/32)*dt*m1+(9/32)*dt*m2,phi,th+(3/32)*dt*l1+(9/32)*dt*l2,tau1,t+(3/8)*dt);

    sig3=sigma(par,u+(1932/2197)*dt*m1-(7200/2197)*dt*m2+(7296/2197)*dt*m3,t+(12/13)*dt);
    
    l4=grns(par,u+(1932/2197)*dt*m1-(7200/2197)*dt*m2+(7296/2197)*dt*m3,phi,th+(1932/2197)*dt*l1-(7200/2197)*dt*l2+(7296/2197)*dt*l3,sig3,tau1,t+(12/13)*dt);
    m4=hrns(par,u+(1932/2197)*dt*m1-(7200/2197)*dt*m2+(7296/2197)*dt*m3,phi,th+(1932/2197)*dt*l1-(7200/2197)*dt*l2+(7296/2197)*dt*l3,tau1,t+(12/13)*dt);

    sig4=sigma(par,u+(439/216)*dt*m1-8*dt*m2+(3680/513)*dt*m3-(845/4104)*dt*m4,t+dt);
    
    l5=grns(par,u+(439/216)*dt*m1-8*dt*m2+(3680/513)*dt*m3-(845/4104)*dt*m4,phi,th+(439/216)*dt*l1-8*dt*l2+(3680/513)*dt*l3-(845/4104)*dt*l4,sig4,tau1,t+dt);
    m5=hrns(par,u+(439/216)*dt*m1-8*dt*m2+(3680/513)*dt*m3-(845/4104)*dt*m4,phi,th+(439/216)*dt*l1-8*dt*l2+(3680/513)*dt*l3-(845/4104)*dt*l4,tau1,t+dt);

    sig5=sigma(par,u-(8/27)*dt*m1+2*dt*m2-(3544/2565)*dt*m3+(1859/4104)*dt*m4-(11/40)*dt*m5,t+0.5*dt);
    
    l6=grns(par,u-(8/27)*dt*m1+2*dt*m2-(3544/2565)*dt*m3+(1859/4104)*dt*m4-(11/40)*dt*m5,phi,th-(8/27)*dt*l1+2*dt*l2-(3544/2565)*dt*l3+(1859/4104)*dt*l4-(11/40)*dt*l5,sig5,tau1,t+0.5*dt);
    m6=hrns(par,u-(8/27)*dt*m1+2*dt*m2-(3544/2565)*dt*m3+(1859/4104)*dt*m4-(11/40)*dt*m5,phi,th-(8/27)*dt*l1+2*dt*l2-(3544/2565)*dt*l3+(1859/4104)*dt*l4-(11/40)*dt*l5,tau1,t+0.5*dt);
    

     
    verror=dt*abs((16/135-25/216)*[l1;m1]+(6656/12825-1408/2565)*[l3;m3]+(28561/56430-2197/4104)*[l4;m4]-(9/50-1/5)*[l5;m5]+(2/55)*[l6;m6]);
    error=max(verror);
    
    if error<=par.err0
        passtep=1;
        dt=min([par.safe*dt*((par.err0/error)^(1/4)) par.dtmax]);
    else
        dt=max([par.safe*dt*((par.err0/error)^(1/4)) par.dtmin]);
    end
       
    
end

if iterrk==par.maxiterrk
    disp([' RK algo did not find solution to sufficient accuracy, error = ',num2str(error)]);
end

th=th+dt*((16/135)*l1+(6656/12825)*l3+(28561/56430)*l4-(9/50)*l5+(2/55)*l6);
u=u+dt*((16/135)*m1+(6656/12825)*m3+(28561/56430)*m4-(9/50)*m5+(2/55)*m6);

sig=sigma(par,u,t+dt);
%disp(l1)
end


function g=grns(par,u,phi,th,sig,tau1,t)
[dbi,dbc]=evoldyke(par,t);
sigm.n=dbc*par.dsigmn;
sigm.t=dbc*par.dsigmt;
%disp('G')
phi=phinew(par,phi,th,u,tau1,t,sigm);
sig=sigma(par,u,t); 
g=exp(-th)-exp(phi)-par.alpha*(par.k*par.b*par.s*(par.c*exp(phi)-1)+par.dsigmn*dbi)/(par.b*(sig));
end

function h=hrns(par,u,phi,th,tau1,t)
[dbi,dbc]=evoldyke(par,t);
sigm.n=dbc*par.dsigmn;
sigm.t=dbc*par.dsigmt;
%disp('H')
phi=phinew(par,phi,th,u,tau1,t,sigm);
h=exp(phi);
end

function sig=sigma(par,u,t) % contrainte normale
[dbi,dbc]=evoldyke(par,t);
sigm.n=dbc*par.dsigmn;
%sigm.t=dbc*par.dsigmt;
sig=1+par.k*(u*par.c-t-par.u0)*par.s+sigm.n; % w = épaisseur cumulée des dikes, à enlever ; ajouter ici terme sig_n
end

function phi=phinew(par,phi,th,u,tau1,t,sigm) % ajouter terme de shear stress (attention conditions)
if par.s==0
    errphi=1;iternew=0;
    while ((errphi>par.err0) & (iternew<par.maxiternew))
        iternew=iternew+1;
        phi0=phi;
        vphi(iternew)=phi;
        f=(par.beta*exp(phi)+(par.f0+th+par.r*phi)*(1+sigm.n))*par.c-par.c/par.b*sigm.t+par.k*(par.c*u-t-par.u0); 
        fp=par.c*(par.r*(1+sigm.n)+par.beta*exp(phi));
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
    M=1+par.b*par.k*(par.c*u-t-par.u0)*par.s;
    %disp(u);
    if M<0
        fm=par.r*M*log(-par.r*M/par.beta)+par.beta^2/(M*par.r)+M*(par.f0+th)-par.s/par.b+(M-1)*par.c/(par.b*par.s)-tau1/par.b-sigm.t*par.c+sigm.n*par.s;
        disp(['error : M < 0 ',num2str(M),' fm : ',num2str(fm)]);
        pause
    else
        while ((errphi>par.err0) & (iternew<par.maxiternew))
            iternew=iternew+1;
            phi0=phi;
            f=par.c*((par.beta*exp(phi)+(M+sigm.n)*(par.f0+th+par.r*phi)))-1/par.b*(sigm.t*par.c-sigm.n*par.s)+(M-1)/(par.b*par.s);
            %disp(f)
            fp=(par.r*(M+sigm.n)+par.beta*exp(phi))*par.c;
            phi=phi-f/fp;
            errphi=abs(phi-phi0);
        end
        if iternew==par.maxiternew
            disp([' NR algo did not find solution to sufficient accuracy, error: ',num2str(errphi),' ',num2str(iternew),' iterations']);
            pause
        end
    end
end
%disp(f)
end

function fco=friction(par,phi,th)
fco=par.b*(par.f0+par.r*phi+th);
end


function [dbi,dbc]=evoldyke(par,t) % t = time dans le loop ; tdike = time du dike ; bmax = rate d'ouverture maximal ; nbd = nombre de jours (en secondes) sur lesquels s'effectue l'ouverture
par.sigd=par.nbd;

dbi=exp(-1/2*((t-(par.tdike+par.delayd))/par.sigd)^2)/(par.sigd*sqrt(2*pi))*par.wd;

dbc=(sqrt(pi/2)*par.sigd*erf((t-(par.tdike+par.delayd))/(sqrt(2)*par.sigd))/(par.sigd*sqrt(2*pi))+1/2)*par.wd;
end

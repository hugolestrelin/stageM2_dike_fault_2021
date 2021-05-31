clear all 
close all

year=3600*24*365;  %--number of seconds in one year
mag_lens_act=1; 
lithos_act=0;
dyke_act=0;
fault_act=0;
code_patin_ressort=1;
P=1e8/(year); % intinial overpressure rate
run sig_defo_v2

%---------------%
%-D parameters--%
%---------------%
param.L=Lf;    %--normal fault length (m)
param.gamma=thetaf-90; %--normal fault angle (deg)
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
param.u0=-1.0;  %--initial displacement (controls initial (compressive) stress)
param.dtau1=0.0*param.rho*param.g*param.z; %--ça je me souviens plus ce que c'est
param.niter=100000; %--number of iterations
param.tfinal=60*year; %--temps de modélisation en année
param.tdike=45*year;
param.nbd=1/365*year;

format long

%--------------------%
%-Termes de forçage--%
%--------------------%

%sigm.n=sigm.n/(param.rho*param.g*param.z);
%sigm.t=sigm.t/(param.rho*param.g*param.z);

sigm0i.n=dsigm.n/(param.rho*param.g*param.z)*param.dc/param.vp;
sigm0i.t=dsigm.t/(param.rho*param.g*param.z)*param.dc/param.vp;

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
sig=1-par.k*par.u0*par.s;%+sigm0i.n*1;  %-Norm. normal stress sig/(rho*g*z)
us=0.0; 
tau1=0.0;
phi=-th;        %-Initial guess for phi (log normalized slip rate)
phi=phinew(par,phi,th,u,tau1,t,sigm0i);  %-compute initial slip rate (from stress balance)

fco=friction(par,phi,th);  %friction coefficient
tau=fco*sig;    %-Norm. shear stress tau/(rho*g*z)

sigm.n=0;%sigm0i.n;
sigm.t=0;%sigm0i.t;
sigmOP=(-(param.rho*param.g*param.z)-(1000*5000*param.g))/(param.rho*param.g*param.z);%1000*5000=Pwater
OPdike=P*year/(param.rho*param.g*param.z)*11.5;%*param.dc/param.vp;
dikecountdown=0;
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

flaguprev=1;
tdikeact=0;
evolsigop=0;
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
    uprev=u;
    phiex=phi;thex=th;uex=u;sigex=sig;tauex=tau;
    [phi,th,u,sig,dt,error]=rk4fehlberg(dt,par,phi,th,u,sig,tau1,t,sigm);
    phi=phinew(par,phi,th,u,tau1,t,sigm);
    fco=friction(par,phi,th);
    tau=fco*sig;
%     if flaguprev>0
%         if uprev>u
%             flaguprev=0;
%             mag_lens_act=1;
%             run sig_defo.m
%         end
%     end
%     
%     if flaguprev==0
%         if uprev>u
%             mag_lens_act=1;
%             fault_act=1;
%             run sig_defo.m
%         end
%     end

    t=t+dt;
    
    if flaguprev>0
        if t>2/3*par.tfinal%uprev<1e5*u
            flaguprev=0;
            disp('premier seisme ?')
        end
    end

%     if tdikeact==0
%         if t>tdike
%             sigm0=Pdike*2;
%             tdikeact=1;
%         end
%     end 

    if flaguprev==0
        
            
        else
            
        end
        if max(sigOP)>OPdike
            %exponential evol -- dike
            %'SIG0 > PDIKE'
            if dikecountdown==0
                dyke_act=1;
                mag_lens_act=1; 
                run sig_defo_v2
                tdike=t;
                dikecountdown=t+3*(par.nbd);
            end
            sigm.n=sigm.n+dsigm.n/(param.rho*param.g*param.z)*b*evoldyke(t,par.nbd,tdike+1/365*year*param.vp/param.dc)*dt*param.dc/param.vp;
            sigm.t=sigm.t+dsigm.t/(param.rho*param.g*param.z)*b*evoldyke(t,par.nbd,tdike+1/365*year*param.vp/param.dc)*dt*param.dc/param.vp;
            if t>dikecountdown
                %retour à la normale
                disp('UN DIKE !')
                dyke_act=0;
                mag_lens_act=1;
                sigm0=P/(param.rho*param.g*param.z)*param.dc/param.vp*dt;
                dikecountdown=0;
                dsigm.n=sigm0i.n;
                dsigm.t=sigm0i.t;
                sigOP=sigOP+dsigOP/(param.rho*param.g*param.z)*b*dt;% +tecto sous dsigop dike slmnt
            end
        %else if FAILLE (cf. dérivée
%             fault_act=1
%             run sig_defo_v2
%             sigOP=sigOP+dsigOP/(param.rho*param.g*param.z)*??*dt;
        else
            %linear evol -- sill
            sigOP=sigOP+dsigOP/(param.rho*param.g*param.z)*param.dc/param.vp*dt;% +tecto sous dsigop sill tecto = par.k*par.u0 ?
            sigm.n=sigm.n+dsigm.n/(param.rho*param.g*param.z)*dt*param.dc/param.vp;
            sigm.t=sigm.t+dsigm.t/(param.rho*param.g*param.z)*dt*param.dc/param.vp; 
        end
        
        
    end
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


function [phi,th,u,sig,dt,error]=rk4fehlberg(dt,par,phi,th,u,sig,tau1,t,sigm)

passtep=0;
iterrk=0;

while ((passtep==0) & (iterrk<=par.maxiterrk))
    
    iterrk=iterrk+1;

    l1=grns(par,u,phi,th,sig,tau1,t,sigm);
    m1=hrns(par,u,phi,th,tau1,t,sigm);
    
    sig1=sigma(par,u+0.25*dt*m1,t+0.25*dt,sigm);
    
    l2=grns(par,u+0.25*dt*m1,phi,th+0.25*dt*l1,sig1,tau1,t+0.25*dt,sigm);
    m2=hrns(par,u+0.25*dt*m1,phi,th+0.25*dt*l1,tau1,t+0.25*dt,sigm);

    sig2=sigma(par,u+(3/32)*dt*m1+(9/32)*dt*m2,t+(3/8)*dt,sigm);
    
    l3=grns(par,u+(3/32)*dt*m1+(9/32)*dt*m2,phi,th+(3/32)*dt*l1+(9/32)*dt*l2,sig2,tau1,t+(3/8)*dt,sigm);
    m3=hrns(par,u+(3/32)*dt*m1+(9/32)*dt*m2,phi,th+(3/32)*dt*l1+(9/32)*dt*l2,tau1,t+(3/8)*dt,sigm);

    sig3=sigma(par,u+(1932/2197)*dt*m1-(7200/2197)*dt*m2+(7296/2197)*dt*m3,t+(12/13)*dt,sigm);
    
    l4=grns(par,u+(1932/2197)*dt*m1-(7200/2197)*dt*m2+(7296/2197)*dt*m3,phi,th+(1932/2197)*dt*l1-(7200/2197)*dt*l2+(7296/2197)*dt*l3,sig3,tau1,t+(12/13)*dt,sigm);
    m4=hrns(par,u+(1932/2197)*dt*m1-(7200/2197)*dt*m2+(7296/2197)*dt*m3,phi,th+(1932/2197)*dt*l1-(7200/2197)*dt*l2+(7296/2197)*dt*l3,tau1,t+(12/13)*dt,sigm);

    sig4=sigma(par,u+(439/216)*dt*m1-8*dt*m2+(3680/513)*dt*m3-(845/4104)*dt*m4,t+dt,sigm);
    
    l5=grns(par,u+(439/216)*dt*m1-8*dt*m2+(3680/513)*dt*m3-(845/4104)*dt*m4,phi,th+(439/216)*dt*l1-8*dt*l2+(3680/513)*dt*l3-(845/4104)*dt*l4,sig4,tau1,t+dt,sigm);
    m5=hrns(par,u+(439/216)*dt*m1-8*dt*m2+(3680/513)*dt*m3-(845/4104)*dt*m4,phi,th+(439/216)*dt*l1-8*dt*l2+(3680/513)*dt*l3-(845/4104)*dt*l4,tau1,t+dt,sigm);

    sig5=sigma(par,u-(8/27)*dt*m1+2*dt*m2-(3544/2565)*dt*m3+(1859/4104)*dt*m4-(11/40)*dt*m5,t+0.5*dt,sigm);
    
    l6=grns(par,u-(8/27)*dt*m1+2*dt*m2-(3544/2565)*dt*m3+(1859/4104)*dt*m4-(11/40)*dt*m5,phi,th-(8/27)*dt*l1+2*dt*l2-(3544/2565)*dt*l3+(1859/4104)*dt*l4-(11/40)*dt*l5,sig5,tau1,t+0.5*dt,sigm);
    m6=hrns(par,u-(8/27)*dt*m1+2*dt*m2-(3544/2565)*dt*m3+(1859/4104)*dt*m4-(11/40)*dt*m5,phi,th-(8/27)*dt*l1+2*dt*l2-(3544/2565)*dt*l3+(1859/4104)*dt*l4-(11/40)*dt*l5,tau1,t+0.5*dt,sigm);
    

     
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

sig=sigma(par,u,t+dt,sigm);

end


function g=grns(par,u,phi,th,sig,tau1,t,sigm)
phi=phinew(par,phi,th,u,tau1,t,sigm);
sig=sigma(par,u,t,sigm); 
g=exp(-th)-exp(phi)-par.alpha*par.k*(par.c*exp(phi)-1)*par.s/sig;
end

function h=hrns(par,u,phi,th,tau1,t,sigm)
phi=phinew(par,phi,th,u,tau1,t,sigm);
h=exp(phi);
end

function sig=sigma(par,u,t,sigm) % contrainte normale
sig=1+par.k*(u*par.c-t-par.u0)*par.s+sigm.n; % w = épaisseur cumulée des dikes, à enlever ; ajouter ici terme sig_n
end

function phi=phinew(par,phi,th,u,tau1,t,sigm) % ajouter terme de shear stress (attention conditions)
if par.s==0
    errphi=1;iternew=0;
    while ((errphi>par.err0) & (iternew<par.maxiternew))
        iternew=iternew+1;
        phi0=phi;
        vphi(iternew)=phi;
        f=par.beta*exp(phi)+par.k*(par.c*u-t-par.u0)-(sigm.t-(par.f0+th+par.r*phi)*(1+sigm.n))*par.c-tau1/par.b; 
        fp=par.r*par.c(1+sigm.n)+par.beta*exp(phi);
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
    if M<0
        fm=par.r*M*log(-par.r*M/par.beta)+par.beta^2/(M*par.r)+M*(par.f0+th)-par.s/par.b+(M-1)*par.c/(par.b*par.s)-tau1/par.b-sigm.t*par.c+sigm.n*par.s;
        disp(['error : M < 0 ',num2str(M),' fm : ',num2str(fm)]);
        pause
    else
        while ((errphi>par.err0) & (iternew<par.maxiternew))
            iternew=iternew+1;
            phi0=phi;
            f=par.beta*exp(phi)+(M+1+sigm.n)*par.c*(par.f0+th+par.r*phi)-par.s/par.b+(M-1)/(par.b*par.s)-sigm.t*par.c+sigm.n*par.s-tau1/par.b;
            fp=par.r*(M+1+sigm.n)*par.c+par.beta*exp(phi);
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


function db=evoldyke(t,nbd,tdike) % t = time dans le loop ; tdike = time du dike ; bmax = rate d'ouverture maximal ; nbd = nombre de jours (en secondes) sur lesquels s'effectue l'ouverture
sigd=4*nbd;
db=exp(-1/2*((t-tdike)/sigd)^2)/(sigd*sqrt(2*pi));
end
clear all 
close all

year=3600*24*365;  %--number of seconds in one year
param.mag_lens_act=1; 
param.dyke_act=0;
param.fault_act=0;
code_patin_ressort=1;
b=0;
P=1e6/(year); % intinial overpressure rate
Pd_ref=P;
run sig_defo_v3
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
param.a=0.001;  %--friction param. 1 (rate-and-state)
param.b=0.002;  %--friction param. 2 (rate-and-state)
param.dc=0.1e-3;%--critical slip (m)
param.vp=0.1/(365*24*3600);%1e-9;  %--tectonic extension rate (m/s)
param.alpha=0.3;%--Linker-Dieterich coeff.
param.u0=0;  %--initial displacement (controls initial (compressive) stress)
param.dtau1=0.0*param.rho*param.g*param.z; %--ça je me souviens plus ce que c'est
param.niter=100000; %--number of iterations
param.tfinal=6*year; %--temps de modélisation en année
param.tdike=1e500*year; % initialisation très (TRÈS) tardive afin de ne pas avoir d'effet sur la loi d'évolution du dike initialement
param.nbd=3*year/120;
param.OPdike=0;%P*year;
param.flaguprev=2*year; % activation de l'intéraction des elmts. Ideal : à automatiser. In pratique : donner une date choisie, attention à la faire varier selon param initiaux.

format long
 
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
%par.delayd=1*year*param.vp/param.dc;
par.wd=b/param.dc;
par.dsigmn_ref=dsigm.n/(param.rho*param.g*param.z)*param.dc/param.vp;% ?
par.dsigmt_ref=dsigm.t/(param.rho*param.g*param.z)*param.dc/param.vp;%/(param.rho*param.g*param.z)?
par.dsigmn=par.dsigmn_ref;
par.dsigmt=par.dsigmt_ref;
par.OPdike=param.OPdike/(param.rho*param.g*param.z);
par.dsigOP_ref=dsigOP/(param.rho*param.g*param.z)*param.dc/param.vp;%/(param.rho*param.g*param.z)
par.dsigOP=par.dsigOP_ref;
par.flaguprev=param.flaguprev*param.vp/param.dc;
par.newsill=par.flaguprev;
par.ndike=0;
par.cumulwidth=0;

disp(['Normalized stiffness : ',num2str(par.k),' Normalized critical stiffness : ',num2str(1-par.r)]);

save parametres.mat param par

%--------------------%
%-Comput parameters--%
%--------------------%
par.err0=1e-10;
par.dtmin=1e-10;
par.dtmax=1e1;
par.maxiterrk=10;
par.maxiternew=30000;
par.safe=0.8;
repres='./results/';

%---------------------%
%-Initial conditions--%
%---------------------%
sigm.n=0;
sigm.t=0;
t=0;            %-Norm. time tvp/dc
dt=1e-6;        %-Norm. Time step
th=log(0.001);  %-Norm log state var
u=0;            %-Norm. slip u/dc
sig=1-par.k*par.b*par.u0*par.s+sigm.n;  %-Norm. normal stress sig/(rho*g*z)
us=0.0; 
tau1=0.0;
phi=-th;        %-Initial guess for phi (log normalized slip rate)
phi=phinew(par,phi,th,u,tau1,t,sigm,param);  %-compute initial slip rate (from stress balance)
dbc=0;
dpc=0;

fco=friction(par,phi,th);  %friction coefficient
tau=fco*sig;    %-Norm. shear stress tau/(rho*g*z)
tpsi=1;
tpsii=[1:10];
dikecountdown=-1;
sigOP_ref=param.rho*param.g*z0/(param.rho*param.g*param.z)*par.dsigOP./par.dsigOP; %z0 profondeur du sill
POP=0;
cumuldike=0;
totaldikewidth=0;
cumulsigop=sigOP_ref;

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
fidfco=fopen([repres,'fco.data'],'w+');
fiddbc=fopen([repres,'dbc.data'],'w+');
fiddpc=fopen([repres,'dpc.data'],'w+');
fidtdw=fopen([repres,'tdw.data'],'w+');
fidsop=fopen([repres,'sop.data'],'w+');

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
    fwrite(fidfco,fco,'real*8');
    fwrite(fiddbc,dbc,'real*8');
    fwrite(fiddpc,dpc,'real*8');
    fwrite(fidtdw,totaldikewidth,'real*8');
    fwrite(fidsop,sigOP,'real*8');
    %----------------------%
    %-Update var.----------%
    %----------------------%
    if param.dyke_act==1
        %disp('on cherche ça ?')
    end
    phiex=phi;thex=th;uex=u;sigex=sig;tauex=tau;
    [phi,th,u,sig,dt,error]=rk4fehlberg(dt,par,phi,th,u,sig,tau1,t,param);
    phi=phinew(par,phi,th,u,tau1,t,sigm,param);
    fco=friction(par,phi,th);
    tau=fco*sig;
    t=t+dt;
    POP=POP+Pd_ref*dt*param.dc/param.vp;
    %par.POP=POP/(param.rho*param.g*param.z);
    
    if t>par.flaguprev
        if max(sigOP)>par.OPdike && param.dyke_act==0 % calculer deux dsigm : dike (pa/m) et decompression sill (pa instantané) car pas m^ evol (donc adimensionalisation)
            disp('Oh le joli dike !')
            param.mag_lens_act=1; 
            P=-POP; %modulation de la decompression en fonction du temps passé
            par.POP=P/(param.rho*param.g*param.z);
            where_dyke=X(1,ixOP(find(sigOP>par.OPdike)),1); % partie où l'on trouve l'endroit où se déclenche dike
            disp(where_dyke)
            run sig_defo_v3
            par.dsigmn_sbh=dsigm.n/(param.rho*param.g*param.z);%*(param.rho*param.g*param.z); % s'annule car dsig / Pa d'OP % est-on sûr de ça ? oui plus ou moins; hmmm probabement que non finalement...
            par.dsigmt_sbh=dsigm.t/(param.rho*param.g*param.z);
            par.dsigOP_sbh=dsigOP/(param.rho*param.g*param.z);
            POP=0;
            param.dyke_act=1;
            param.mag_lens_act=0; 
            b=0.01*par.POP/(-0.212609338316523); % wd proportionnel à POP sur la base de la première valeur de POP ; est-ce juste ?
            par.wd=b/param.dc;
            run sig_defo_v3
            par.tdike=t;
            dikecountdown=t+3*3*(par.nbd);
            par.ndike=par.ndike+1;
            par.dsigmn_dike=dsigm.n/(param.rho*param.g*param.z)*param.dc;
            par.dsigmt_dike=dsigm.t/(param.rho*param.g*param.z)*param.dc;
            par.dsigOP_dike=dsigOP/(param.rho*param.g*param.z)*param.dc;
            par.newsill=dikecountdown;
        else
            if t<dikecountdown
                param.dyke_act=1;
            else 
                if param.dyke_act==1 %&& par.ndike>1
                    par.cumulwidth=par.cumulwidth+par.wd;
                end
                %retour à la normale
                param.dyke_act=0;
                param.mag_lens_act=1;
                %sigm0=P/(param.rho*param.g*param.z)*param.dc/param.vp*dt;
                %dikecountdown=0;
                dsigm.n=par.dsigmn_ref;
                dsigm.t=par.dsigmt_ref;
                par.dsigmn=dsigm.n;
                par.dsigmt=dsigm.t;
                par.dsigOP=par.dsigOP_ref;
                %sigOP=sigOP+par.dsigOP*b*dt;
            end
%         elseif (u-uprev)/(par.c*dt)>0.015
%             fault_act=1;
%             dyke_act=0;
%             mag_lens_act=0;
%             run sig_defo_v3
%             sigOP=sigOP+par.dsigOP*??*dt;
%         else
%             %linear evol -- sill
%             % +tecto sous dsigop sill tecto = par.k*par.u0 ?
%             sigm.n=sigm.n+dsigm.n/(param.rho*param.g*param.z)*dt*param.dc/param.vp;
%             sigm.t=sigm.t+dsigm.t/(param.rho*param.g*param.z)*dt*param.dc/param.vp; 
%             sigOP=sigOP+par.dsigOP*b*dt;
        end
        [dbi,dbc,dpi,dpc]=evoldyke(param,par,t);
        totaldikewidth=dbc+par.cumulwidth;
        if param.dyke_act==1
            %disp('on est là ?')
            sigm.n=dbc*par.dsigmn_dike+dpc*par.dsigmn_sbh;
            sigm.t=dbc*par.dsigmt_dike+dpc*par.dsigmt_sbh;
            sigOP=dbc*par.dsigOP_dike+dpc*par.dsigOP_sbh+sigOP_ref; %/ simple addition constante des résultats du champs des contraintes non ?
        else
            sigm.n=dpc*par.dsigmn;
            sigm.t=dpc*par.dsigmt;
            sigOP=dpc*par.dsigOP+sigOP_ref;
        end
    end
    
    %sigOP=sigOP+par.dsigOP*dt;
    if t>tpsii(tpsi)/10*par.tfinal
            disp([num2str(floor(100*t/par.tfinal)),' % completed']);
            tpsi=tpsi+1;
            %disp(sigOP)
    end
    %disp(sigOP)
    cumulsigop=[cumulsigop sigOP];
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
fclose(fidfco);
fclose(fiddbc);
fclose(fiddpc);
fclose(fidtdw);
fclose(fidsop);

function [phi,th,u,sig,dt,error]=rk4fehlberg(dt,par,phi,th,u,sig,tau1,t,param)

passtep=0;
iterrk=0;

while ((passtep==0) & (iterrk<=par.maxiterrk))
    
    iterrk=iterrk+1;

    l1=grns(par,u,phi,th,tau1,t,param);
    m1=hrns(par,u,phi,th,tau1,t,param);
        
    l2=grns(par,u+0.25*dt*m1,phi,th+0.25*dt*l1,tau1,t+0.25*dt,param);
    m2=hrns(par,u+0.25*dt*m1,phi,th+0.25*dt*l1,tau1,t+0.25*dt,param);
    
    l3=grns(par,u+(3/32)*dt*m1+(9/32)*dt*m2,phi,th+(3/32)*dt*l1+(9/32)*dt*l2,tau1,t+(3/8)*dt,param);
    m3=hrns(par,u+(3/32)*dt*m1+(9/32)*dt*m2,phi,th+(3/32)*dt*l1+(9/32)*dt*l2,tau1,t+(3/8)*dt,param);
    
    l4=grns(par,u+(1932/2197)*dt*m1-(7200/2197)*dt*m2+(7296/2197)*dt*m3,phi,th+(1932/2197)*dt*l1-(7200/2197)*dt*l2+(7296/2197)*dt*l3,tau1,t+(12/13)*dt,param);
    m4=hrns(par,u+(1932/2197)*dt*m1-(7200/2197)*dt*m2+(7296/2197)*dt*m3,phi,th+(1932/2197)*dt*l1-(7200/2197)*dt*l2+(7296/2197)*dt*l3,tau1,t+(12/13)*dt,param);
    
    l5=grns(par,u+(439/216)*dt*m1-8*dt*m2+(3680/513)*dt*m3-(845/4104)*dt*m4,phi,th+(439/216)*dt*l1-8*dt*l2+(3680/513)*dt*l3-(845/4104)*dt*l4,tau1,t+dt,param);
    m5=hrns(par,u+(439/216)*dt*m1-8*dt*m2+(3680/513)*dt*m3-(845/4104)*dt*m4,phi,th+(439/216)*dt*l1-8*dt*l2+(3680/513)*dt*l3-(845/4104)*dt*l4,tau1,t+dt,param);
    
    l6=grns(par,u-(8/27)*dt*m1+2*dt*m2-(3544/2565)*dt*m3+(1859/4104)*dt*m4-(11/40)*dt*m5,phi,th-(8/27)*dt*l1+2*dt*l2-(3544/2565)*dt*l3+(1859/4104)*dt*l4-(11/40)*dt*l5,tau1,t+0.5*dt,param);
    m6=hrns(par,u-(8/27)*dt*m1+2*dt*m2-(3544/2565)*dt*m3+(1859/4104)*dt*m4-(11/40)*dt*m5,phi,th-(8/27)*dt*l1+2*dt*l2-(3544/2565)*dt*l3+(1859/4104)*dt*l4-(11/40)*dt*l5,tau1,t+0.5*dt,param);
     
    verror=dt*abs((16/135-25/216)*[l1;m1]+(6656/12825-1408/2565)*[l3;m3]+(28561/56430-2197/4104)*[l4;m4]-(9/50-1/5)*[l5;m5]+(2/55)*[l6;m6]);
    error=max(verror);
    
    if error<=par.err0
        passtep=1;
        dt=min([par.safe*dt*((par.err0/error)^(1/4)) par.dtmax]);
        %disp('gros dt ?')
        %disp(dt)
        %disp(error)
    else
        dt=max([par.safe*dt*((par.err0/error)^(1/4)) par.dtmin]);
        %disp('petit dt ?')
        %disp(dt)
    end
       
    
end

if iterrk==par.maxiterrk
    disp([' RK algo did not find solution to sufficient accuracy, error = ',num2str(error)]);
end

th=th+dt*((16/135)*l1+(6656/12825)*l3+(28561/56430)*l4-(9/50)*l5+(2/55)*l6);
u=u+dt*((16/135)*m1+(6656/12825)*m3+(28561/56430)*m4-(9/50)*m5+(2/55)*m6);

sig=sigma(par,u,t+dt,param);
%disp(l1)
end


function g=grns(par,u,phi,th,tau1,t,param)
[dbi,dbc,dpi,dpc]=evoldyke(param,par,t);
    if param.dyke_act==1
        sigm.n=dbc*par.dsigmn_dike+dpc*par.dsigmn_sbh;
        sigm.t=dbc*par.dsigmt_dike+dpc*par.dsigmt_sbh;
    else
        sigm.n=dpc*par.dsigmn;
        sigm.t=dpc*par.dsigmt;
    end
phi=phinew(par,phi,th,u,tau1,t,sigm,param);
sig=sigma(par,u,t,param); 
    if param.dyke_act==1
        g=exp(-th)-exp(phi)-par.alpha*(par.k*par.b*par.s*(par.c*exp(phi)-1+dbi)+par.dsigmn_dike*dbi+par.dsigmn_sbh*dpi)/(par.b*(sig)); % ?? NECESSITÉ des conditions ? Besoin de l'addition de dpi/dpc ? same Q? pour toutes les autres fct / +dpi ?
    else
        g=exp(-th)-exp(phi)-par.alpha*(par.k*par.b*par.s*(par.c*exp(phi)-1)+(par.dsigmn)*dpi)/(par.b*(sig)); % Ou alors, non sgnificativité de dbc après un certain temps donc addition toujours faisable ? / +dpi ?
        % addition de dpi ? le faire moduler (*par.POP) ? l'enlever ?
        % significtion = ajout d'un incrément de volume du sill ?
        % négligeable car augm pression vient du fait qu'on ne "grossit" ap ?
    end
end

function h=hrns(par,u,phi,th,tau1,t,param)
[dbi,dbc,dpi,dpc]=evoldyke(param,par,t);
    if param.dyke_act==1
        sigm.n=dbc*par.dsigmn_dike+dpc*par.dsigmn_sbh;
        sigm.t=dbc*par.dsigmt_dike+dpc*par.dsigmt_sbh;
    else
        sigm.n=dpc*par.dsigmn;
        sigm.t=dpc*par.dsigmt;
    end
phi=phinew(par,phi,th,u,tau1,t,sigm,param);
h=exp(phi);
end

function sig=sigma(par,u,t,param) % contrainte normale
    [dbi,dbc,dpi,dpc]=evoldyke(param,par,t);
    if param.dyke_act==1
        sigm.n=dbc*par.dsigmn_dike+dpc*par.dsigmn_sbh;
        sigm.t=dbc*par.dsigmt_dike+dpc*par.dsigmt_sbh;
    else
        sigm.n=dpc*par.dsigmn;
        sigm.t=dpc*par.dsigmt;
    end
if param.dyke_act==1
    sig=1+par.k*par.b*(u*par.c-t-par.u0+max(0,par.cumulwidth)+dbc)*par.s+sigm.n; %+dpc ?
else
    sig=1+par.k*par.b*(u*par.c-t-par.u0+max(0,par.cumulwidth))*par.s+sigm.n;
end
%disp(sig)
end

function phi=phinew(par,phi,th,u,tau1,t,sigm,param) % ajouter terme de shear stress (attention conditions)
[dbi,dbc,dpi,dpc]=evoldyke(param,par,t);
if par.s==0
    errphi=1;iternew=0;
    while ((errphi>par.err0) & (iternew<par.maxiternew))
        iternew=iternew+1;
        phi0=phi;
        vphi(iternew)=phi;
        if param.dyke_act==1
            f=(par.beta*exp(phi)+(par.f0+th+par.r*phi)*(1+sigm.n))*par.c-par.c/par.b*sigm.t+par.k*(par.c*u-t-par.u0+max(0,par.cumulwidth)+dbc); %+dpc ?
        else
            f=(par.beta*exp(phi)+(par.f0+th+par.r*phi)*(1+sigm.n))*par.c-par.c/par.b*sigm.t+par.k*(par.c*u-t-par.u0+max(0,par.cumulwidth));
        end
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
    if param.dyke_act==1
        M=1+par.b*par.k*(par.c*u-t-par.u0+max(0,par.cumulwidth)+dbc)*par.s; % +dpc ?
    else
        M=1+par.b*par.k*(par.c*u-t-par.u0+max(0,par.cumulwidth))*par.s;
    end
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
    fco=(par.f0+par.r*phi+th);
end

function [dbi,dbc,dpi,dpc]=evoldyke(param,par,t)
    par.sigd=par.nbd;
    par.delayd=3*par.sigd;
    if param.dyke_act==1
        dbi=exp(-1/2*((t-(par.tdike+par.delayd))/par.sigd)^2)/(par.sigd*sqrt(2*pi))*par.wd;
        dbc=(sqrt(pi/2)*par.sigd*erf((t-(par.tdike+par.delayd))/(sqrt(2)*par.sigd))/(par.sigd*sqrt(2*pi))+1/2)*par.wd;
        dpi=exp(-1/2*((t-(par.tdike+par.delayd))/par.sigd)^2)/(par.sigd*sqrt(2*pi))*par.POP;
        dpc=(sqrt(pi/2)*par.sigd*erf((t-(par.tdike+par.delayd))/(sqrt(2)*par.sigd))/(par.sigd*sqrt(2*pi))+1/2)*par.POP;
    else
        dbc=0;
        dbi=0;
        if t<par.flaguprev
           dpi=0;
           dpc=0;
        else
           dpi=(t-par.newsill);
           dpc=((t^2)/2-par.newsill*t+par.newsill^2/2);
        end
    end
end


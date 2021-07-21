clear all 
close all

code_patin_ressort=1;
year=3600*24*365;  %--number of seconds in one year

%---------------%
%-D parameters--%
%---------------%
param.L=1500/sind(60);    %--normal fault length (m)
param.gamma=60; %--normal fault angle (deg)
param.widthsill=502; %--epaisseur du sill
param.thsill=0; %--orientation du sill /x1
param.dfs=5000; %--distance faille-sill au niveau de l'extrémité la plus profonde de la faille
param.zsill=-1500; %--profondeur du sill, 0 correspond à la surface de l'océan
param.xsill=0; %--position centre sill
param.dd=0; %--profondeur du top du dike
param.Ld=param.zsill; %--longuer du dike
param.z=abs(param.zsill+param.zsill+(param.L*sind(param.gamma)))/2;    %--depth (m) abs((xsiup(2)+xsidown(2))/2)
param.rho=3e3;  %--rock density (kg/m3)
param.g=9.81;   %--gravity (m2/s)
param.mu=3e10;  %--shear modulus (Pa)
param.nu=0.254; %--Poisson coeff
param.f0=0.6;   %--static friction coeff.
param.a=0.001;  %--friction param. 1 (rate-and-state)
param.b=0.002;  %--friction param. 2 (rate-and-state)
param.dc=0.1e-3;%--critical slip (m)
param.vp=0.1/(365*24*3600);%1e-9;  %--tectonic extension rate (m/s)
param.alpha=0.3;%--Linker-Dieterich coeff.
param.u0=0;  %--initial displacement (controls initial (compressive) stress)
param.dtau1=0.0*param.rho*param.g*param.z; %--ça je me souviens plus ce que c'est
param.niter=100000; %--number of iterations
param.pasgrille=10; %--pas de la grille pour le calcul des déformations
param.tfinal=15*year; %--temps de modélisation en année
param.tdike=1e500*year; % initialisation très (TRÈS) tardive afin de ne pas avoir d'effet sur la loi d'évolution du dike initialement
param.nbd=1*year/365;
param.OPdike=0;%P*year;
param.flaguprev=1.15*year; % activation de l'intéraction des elmts. Ideal : à automatiser. In pratique : donner une date choisie, attention à la faire varier selon param initiaux.

format long

%-----------%
% runs de ref
%-----------%
param.dyke_act=0;
param.mag_lens_act=1;
param.fault_act=0;
b=0;
P=-1;
run sig_defo_v3
par.dsigOP_sbh_ref=dsigOP;
par.dsigm_sbh_refn=dsigm.n;%/(param.rho*param.g*param.z);
par.dsigm_sbh_reft=dsigm.t;%/(param.rho*param.g*param.z);

param.dyke_act=1;
param.mag_lens_act=0; 
param.fault_act=0;
P=0;
b=1;
xOP=-param.widthsill/2-3*param.pasgrille:param.pasgrille:param.widthsill/2+3*param.pasgrille;
where_dyke=xOP(find(xOP==param.xsill-param.widthsill/2-2*param.pasgrille));
run sig_defo_v3
par.dsigmn_dike_ref=dsigm.n/(param.rho*param.g*param.z)*param.dc;
par.dsigmt_dike_ref=dsigm.t/(param.rho*param.g*param.z)*param.dc;
par.dsigOP_dike_ref=dsigOP/(param.rho*param.g*param.z)*param.dc;

param.mag_lens_act=0; 
param.dyke_act=0;
param.fault_act=1;
run sig_defo_v3
% BANCAL ?
par.dsigOPfault=(dsigOP/(param.rho*param.g*param.z))*param.dc;

param.mag_lens_act=1; 
param.dyke_act=0;
param.fault_act=0;
b=0;
P=0.4;%1e7/(year); % intinial overpressure rate
Pd_ref=P;
run sig_defo_v3

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
tpsi2=1;
tpsii=[1:10];
isave=1000;
tpsii2=[1:isave];
dikecountdown=-1;
sigOP_ref=param.rho*param.g*param.zsill/(param.rho*param.g*param.z)*ones(size(par.dsigOP)); %z0 profondeur du sill
sigOP=sigOP_ref;
POP=0;
cumuldike=0;
totaldikewidth=0;
cumulsigop=sigOP_ref;
dike_deficit=0;

catdike.width=[];
catdike.t=[];
catdike.n=[];
catdike.nbd=[];
catdike.where=[];

catseisme.t=[];

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
fidpop=fopen([repres,'pop.data'],'w+');

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
    fwrite(fidpop,POP,'real*8');
    %----------------------%
    %-Update var.----------%
    %----------------------%
    phiex=phi;thex=th;uex=u;sigex=sig;tauex=tau;
    [phi,th,u,sig,dt,error]=rk4fehlberg(dt,par,phi,th,u,sig,tau1,t,param);
    phi=phinew(par,phi,th,u,tau1,t,sigm,param);
    fco=friction(par,phi,th);
    tau=fco*sig;
    %par.POP=POP/(param.rho*param.g*param.z);
    if exp(phi)>1e5
        catseisme.t=[catseisme.t t];
    end
    if t>par.flaguprev
        dike_deficit=dike_deficit+par.k*par.b*dt;
        %if max(sigOP)>par.OPdike && param.dyke_act==0 % calculer deux dsigm : dike (pa/m) et decompression sill (pa instantané) car pas m^ evol (donc adimensionalisation)
        if sigOP(find(xOP==param.xsill-param.widthsill/2-2*param.pasgrille))>par.OPdike && param.dyke_act==0
            disp('Oh le joli dike !')
            par.POP=-POP/(param.rho*param.g*param.z);
            where_dyke=xOP(find(xOP==param.xsill-param.widthsill/2-2*param.pasgrille));
            par.dsigmn_sbh=par.dsigm_sbh_refn;%*(param.rho*param.g*param.z); % s'annule car dsig / Pa d'OP % est-on sûr de ça ? oui plus ou moins; hmmm probabement que non finalement...
            par.dsigmt_sbh=par.dsigm_sbh_reft;
            par.dsigOP_sbh=par.dsigOP_sbh_ref;
            POP=0;
            % BANCAL
            %b=abs(par.POP*(param.rho*param.g*param.z))*param.widthsill^2/abs(param.zsill)*(1-param.nu)/param.mu;%0.01*par.POP/(-0.212609338316523); % wd proportionnel à POP sur la base de la première valeur de POP ; est-ce juste ?
            %b=dike_deficit*param.dc; % première approx NR
            %disp(b)
            %SIG_dyke_NR=dyke(param.dd,param.Ld,0,param.zsill+param.pasgrille,b,param.nu,param.mu);
            SIG_dyke_NR=dyke(param.dd,param.Ld,0,param.zsill+param.pasgrille,1,param.nu,param.mu);
            b=dike_deficit*(param.rho*param.g*param.z)/(abs(SIG_dyke_NR.xx));
            %NR_cond=(dike_deficit-abs(SIG_dyke_NR.xx)/(param.rho*param.g*param.z)*param.dc);
            %disp(NR_cond)
%             if NR_cond>1e-9
%                 while NR_cond>0
%                     db=1e-8*b;
%                     db_SIG_dyke_NR=dyke(param.dd,param.Ld,0,param.zsill+param.pasgrille,b+db,param.nu,param.mu);
%                     dSIG_dyke_NR=(db_SIG_dyke_NR.xx-SIG_dyke_NR.xx)/(param.rho*param.g*param.z*db)*param.dc; %db/param.dc ??
% %                     if SIG_dyke_NR.xx/(param.rho*param.g*param.z)*param.dc/dSIG_dyke_NR>0.02
% %                         b=b+0.01;
% %                     else
% %                         b=b+SIG_dyke_NR.xx/(param.rho*param.g*param.z)*param.dc/dSIG_dyke_NR;  % - * - = + ?
% %                     end
%                     bprev=b;
%                     b=b+SIG_dyke_NR.xx/(param.rho*param.g*param.z)*param.dc/dSIG_dyke_NR;  % - * - = + ?
%                     SIG_dyke_NR=dyke(param.dd,param.Ld,0,param.zsill+param.pasgrille,b,param.nu,param.mu);
%                     %NR_cond=(sigOP(find(sigOP==max(sigOP)))-SIG_dyke_NR.xx/(param.rho*param.g*param.z)*param.dc);
%                     NR_cond=(dike_deficit-abs(SIG_dyke_NR.xx)/(param.rho*param.g*param.z)*param.dc);
%                     %disp(NR_cond)
%                     %disp(b)
%                     %disp('----')
%                 end
%             end
            %b=bprev;%-SIG_dyke_NR.xx/(param.rho*param.g*param.z)*param.dc/dSIG_dyke_NR;
            %disp(bprev)
            param.dyke_act=1;
            disp(b)
            par.wd=b/param.dc;
            par.tdike=t;
            % BANCAL
            dikecountdown=t+5.6*(par.nbd);
            dike_deficit_prev=dike_deficit;
            dike_deficit=0;
            par.ndike=par.ndike+1;
            par.dsigmn_dike=par.dsigmn_dike_ref;
            par.dsigmt_dike=par.dsigmt_dike_ref;
            par.dsigOP_dike=par.dsigOP_dike_ref;
            %disp(dsigOP)
            par.newsill=dikecountdown;
            catdike.width=[catdike.width par.wd];
            catdike.t=[catdike.t t];
            %catdike.nbd=[catdike.nbd par.nbd];
            catdike.n=[catdike.n par.ndike];
            %catdike.where=[catdike.where where_dyke];
        else
            if t<dikecountdown
                param.dyke_act=1;
            else 
                if param.dyke_act==1 %&& par.ndike>1
                    par.cumulwidth=par.cumulwidth+par.wd;
                end
                %retour à la normale
                POP=POP+Pd_ref*dt*param.dc/param.vp;
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
        end
        [dbi,dbc,dpi,dpc]=evoldyke(param,par,t);
        totaldikewidth=dbc+par.cumulwidth;
        if param.dyke_act==1
            %disp('on est là ?')
            sigm.n=dbc*par.dsigmn_dike+dpc*par.dsigmn_sbh;
            sigm.t=dbc*par.dsigmt_dike+dpc*par.dsigmt_sbh;
            %sigOP=dbc*par.dsigOP_dike+dpc*par.dsigOP_sbh+sigOP_ref+par.dsigOPfault*exp(phi)*dt; %/ simple addition constante des résultats du champs des contraintes non 
            %sigOP=sigOP+(par.k*par.b+1e2*dbi*par.dsigOP_dike-dpi*par.dsigOP_sbh+par.dsigOPfault*exp(phi))*dt; %/ simple addition constante des résultats du champs des contraintes non ?
            
            %sigOP=-2*ones(size(par.dsigOP));%sigOP+(par.k*par.b-dpi*dike_deficit_prev/par.POP-dpi*par.dsigOP_sbh+par.dsigOPfault*exp(phi))*dt; %/ simple addition constante des résultats du champs des contraintes non ?
            sigOP=sigOP+(par.k*par.b-dpi*dike_deficit_prev/par.POP-dpi*par.dsigOP_sbh+par.dsigOPfault*exp(phi))*dt;
            %disp(sigOP)
            %disp(dbi)
            %disp(dpi)
        else
            sigm.n=dpc*par.dsigmn;
            sigm.t=dpc*par.dsigmt;
            %sigOP=dpc*par.dsigOP+sigOP_ref+par.dsigOPfault*exp(phi)*dt;%*u
            % BANCAL *par.b
            sigOP=sigOP+(par.k*par.b+dpi*par.dsigOP+par.dsigOPfault*exp(phi))*dt;%*u
        end
        %disp(sigOP)
    end
    t=t+dt;
    %sigOP=sigOP+par.dsigOP*dt;
    if t>tpsii(tpsi)/10*par.tfinal
            disp([num2str(floor(100*t/par.tfinal)),' % completed']);
            tpsi=tpsi+1;
            %disp(sigOP)
    end
    if t>=tpsii2(tpsi2)/isave*par.tfinal
            save(['results/' num2str(tpsi2)],'t','sigOP');
            tpsi2=tpsi2+1;
            %cumulsigop=[cumulsigop sigOP];
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
fclose(fidfco);
fclose(fiddbc);
fclose(fiddpc);
fclose(fidtdw);
fclose(fidsop);
fclose(fidpop);

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
           dpi=1;
           dpc=(t-par.newsill);
        end
    end
end


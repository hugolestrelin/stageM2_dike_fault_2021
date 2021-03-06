%--------------%
%--parameters--%
%--------------%
%code_patin_ressort=0;

year=3600*24*365;  %--number of seconds in one year
lithos_act=0;
    
if code_patin_ressort==1
   %disp('interaction activated')
   mag_lens_act=param.mag_lens_act; 
   dyke_act=param.dyke_act;
   fault_act=param.fault_act;
else
    disp('!!! ATTENTION NO INTERACTION AVEC PATIN-RESSORT !!!')
    mag_lens_act=1; 
    dyke_act=0;
    fault_act=0;
end

%--INTERACTIONS PARAMETERS--%

%--SILL--%
if code_patin_ressort==1
    if param.mag_lens_act==1
        P=P;
    else
        P=0;
    end
   L=param.widthsill;
   theta=param.thsill;
   z0=param.zsill;
   dfs=param.dfs;
   g=param.g;
   rho=param.rho;
   pas=param.pasgrille;
   nus=param.nu;
   x0=param.xsill;
else
   P=1e8/(year*1e1);    % initial overpressure rate in the crack (Pa/s)
   L=502;    % crack width (m) ; take an even nbr
   theta=0; % crack angle relative to the x axis (deg) non implement√©
   z0=-1000;   % crack depth (m)
   x0=0; % crack horizontal location (m)
   nus=0.254;   % Poisson coeff
   dfs=1000;%[0,1000,2000,3000,4000,5000,10000,15000,20000]; % distance btw end of the fault and sill (m)
   %--LITHOSTATIC PRESSURE--%
   g=9.8; %--free-fall cte (m/s)
   rho=2900; %--rock density (kg/m3)
   %--BOX--%
   pas=10;
end

%--DIKE--%
if code_patin_ressort==1
    if param.dyke_act==1
        b=b;
    else
        b=0;
    end
    mu=param.mu;
    nud=nus;
    dd=param.dd;
    Ld=param.Ld;
else
    b=0.1; % crack opening (m)
    where_dyke=0;
    mu=3e10; % shear modulus (Pa)
    nud=nus; % poisson coefficient
    dd=0; % depth to top of crack (m) %-100
    Ld=z0; % height of crack (m)
end


sign_tt=dfs;
tau_tt=dfs;
sign_ttt=dfs;
tau_ttt=dfs;
sign_ttt_wd=dfs;
tau_ttt_wd=dfs;
rr=size(dfs);
SIG_lithos={};
SIG_dyke={};
SIG_fault={};
SIG_Mag={};
SIGT={};
sigOP=[];

%--FAULT--%
if code_patin_ressort==1
    thetaf=param.gamma;
    nu=nud;
    Lf=param.L;
else
    thetaf=60;%+180? % fault dip
    % if thetaf==0
    %     Lf=abs(z0);
    % else
    %     Lf=abs(z0)/sind(thetaf); % surface of displacement on fault (?)
    % end
    Lf=500;
    nu=nud;
end

pasxsi=pas;
xsilength=100;

for k=1:rr(2)
    %--BOX--%
    d_fs=dfs(k);
    box=[dfs(k)+15000,abs(z0)+10000];
    x =-box(1):pas:box(1);
    z = -box(2):pas:0;
    [X,Z] = meshgrid(x,z);
    s=size(X);
    tau_tot=X;
    sign_tot=X;
   
    depthf=z0;
    topf=z0+(Lf*sind(thetaf));
    xsidown=[x0+L/2+d_fs,depthf];
    xsiup=[x0+L/2+d_fs+((Lf)*cosd(thetaf)),topf];% origin dislocation %z0+(Lf*sind(thetaf))/2 % 
    df=xsiup(2); % fault top depth +(Lf*sind(thetaf))/2
    if cosd(thetaf)==0
        xfp=zeros([1 10]);
        %xsi1=zeros([1 pasxsi]);
        zfp=depthf:(Lf*sind(thetaf))/9:topf;
        %xsi2=xsi(2)-(xsilength*sind(thetaf))/2:(xsilength*sind(thetaf))/(pasxsi-1):xsi(2)+(xsilength*sind(thetaf))/2;
    elseif sind(thetaf)==0
        xfp=xsidown(1):(Lf*cosd(thetaf))/9:xsiup(1); %!!! SEULEMENT SI FAILLE √Ä DROITE SILL, AJOUTER CONDITION SI NON !!!
        %xsi1=xsi(1)-(xsilength*cosd(thetaf))/2:(xsilength*cosd(thetaf))/(pasxsi-1):xsi(1)+(xsilength*cosd(thetaf))/2;
        zfp=zeros([1 10]);
        %xsi2=zeros([1 pasxsi]);
    else
        xfp=xsidown(1):(Lf*cosd(thetaf))/10:xsiup(1);
        %xsi1=xsi(1)-(xsilength*cosd(thetaf))/2:(xsilength*cosd(thetaf))/(pasxsi-1):xsi(1)+(xsilength*cosd(thetaf))/2;
        zfp=xsidown(2):(Lf*sind(thetaf))/10:xsiup(2);
        %xsi2=xsi(2)-(xsilength*sind(thetaf))/2:(xsilength*sind(thetaf))/(pasxsi-1):xsi(2)+(xsilength*sind(thetaf))/2;
    end
    % where is le milieu de la faille ?
    F=X>abs(xsidown(1)+xsiup(1))/2;
    Arrr=X(F);
    idxx=find(X(1,:,1)==Arrr(1));
    F=Z>(xsidown(2)+xsiup(2))/2;
    Arrr=Z(F);
    idxz=find(Z(:,1,1)==Arrr(1));
    % where is la moyenne de la faille ?
    xfpi=size(xfp);
    idxfp=xfp;
    idzfp=zfp;
    for l=1:xfpi(2)
        idxfp(l)=find(X(1,:,1)>xfp(l)-pas/2 & X(1,:,1)<xfp(l)+pas/2);
        idzfp(l)=find(Z(:,1,1)>zfp(l)-pas/2 & Z(:,1,1)<zfp(l)+pas/2);
    end
    %d√©finition de la ligne d'√©valuation de SIGOP pour savoir quand on peut
    %faire un dike :
    xOP=-L/2-3*pas:pas:L/2+3*pas;%-L/2-2*pas:pas:L/2+2*pas; % demi-sill : 0:pas:L/2+3*pas; %% SI BESOIN EFFECTUER CHANGEMENT TAILLE VECTEUR D'√ČVALUATION DES SIG AU DESSUS DU SILL
    sxop=size(xOP);
    zOP=(z0+pas)*ones(1,sxop(2));
    izOP=find(Z(:,1,1)==zOP(1))*ones(1,sxop(2));
    ixOP=[];
    for l=1:sxop(2)
        ixOP=[ixOP,find(X(1,:,1)>xOP(l)-pas/2 & X(1,:,1)<xOP(l)+pas/2)];
    end
    if lithos_act==1   
        SIG_lithos=lithos(rho,g,z);
    else
        SIG_lithos.xx = zeros(size(X));
        SIG_lithos.zz = zeros(size(X));
        SIG_lithos.xz = zeros(size(X));
    end
    if mag_lens_act==1
        SIG_Mag = mag_lens_f(X,Z,P,L,nus,z0,x0);
    else
        SIG_Mag.xx = zeros(size(X));
        SIG_Mag.zz = zeros(size(X));
        SIG_Mag.xz = zeros(size(X));
    end
    if dyke_act==1
        SIG_dyke = dyke(dd,Ld,X-where_dyke(1),Z,b,nud,mu);
    else
        SIG_dyke.xx = zeros(size(X));
        SIG_dyke.zz = zeros(size(X));
        SIG_dyke.xz = zeros(size(X));
    end
    if fault_act==1
        [SIG_faultup,Uu1,Uu2] = fault_okada_v1(X,Z,mu,Lf,df,nus,thetaf,xsiup);
        [SIG_faultbot,Ub1,Ub2] = fault_okada_v1(X,Z,mu,Lf,df,nus,thetaf,xsidown);
        %SIG_fault=SIG_faultup-SIG_faultbot;
%         SIG_fault.xx = SIG_faultbot.xx-SIG_faultup.xx; 
%         SIG_fault.zz = SIG_faultbot.zz-SIG_faultup.zz;
%         SIG_fault.xz = SIG_faultbot.xz-SIG_faultup.xz;
        SIG_fault.xx = SIG_faultup.xx-SIG_faultbot.xx; 
        SIG_fault.zz = SIG_faultup.zz-SIG_faultbot.zz;
        SIG_fault.xz = SIG_faultup.xz-SIG_faultbot.xz;
    else
        SIG_fault.xx = zeros(size(X));
        SIG_fault.zz = zeros(size(X));
        SIG_fault.xz = zeros(size(X));
        U1=0;
        U2=0;
    end
    SIGT.xx = SIG_Mag.xx+SIG_lithos.xx/abs(-rho*g)+SIG_dyke.xx+SIG_fault.xx;%/(mu*abs(Lf)/(2*pi*(1-nus)))
    SIGT.zz = SIG_Mag.zz+SIG_lithos.zz/abs(-rho*g)+SIG_dyke.zz+SIG_fault.zz;
    SIGT.xz = SIG_Mag.xz+SIG_lithos.xz/abs(-rho*g)+SIG_dyke.xz+SIG_fault.xz;
    tau_tot(:,:)=(SIGT(1).xx-SIGT(1).zz)*cosd(thetaf)*sind(thetaf)+SIGT(1).xz*(sind(thetaf)^2-cosd(thetaf)^2); % rep√®re de pierre
    sign_tot(:,:)=SIGT(1).xx*sind(thetaf)^2+SIGT(1).zz*cosd(thetaf)^2-2*SIGT(1).xz*cosd(thetaf)*sind(thetaf); % rep√®re de pierre
    %tau_tot_wd(:,:)=(SIGT(1).xx-SIGT(1).zz)*cosd(90)*sind(90)+SIGT(1).xz*(sind(90)^2-cosd(90)^2); 
    %sign_tot_wd(:,:)=SIGT(1).xx*sind(90)^2+SIGT(1).zz*cosd(90)^2-2*SIGT(1).xz*cosd(90)*sind(90);
      
    sign_ttt(k)=sign_tot(idxz,idxx);
    tau_ttt(k)=tau_tot(idxz,idxx);
    sign_ttt_wd(k)=sign_tot(idxz,round(s(2)/2));
    %tau_ttt_wd(k)=tau_tot(idxz,round(s(2)/2));
    sigOP=(SIGT.xx(izOP(1),ixOP)+SIGT.zz(izOP(1),ixOP))/2-(SIGT.xz(izOP(1),ixOP));
    sss=0;
    ttt=0;
    for m=2:xfpi(2)
        sss=sss+sign_tot(idzfp(m),idxfp(m));
        ttt=ttt+tau_tot(idzfp(m),idxfp(m));
    end
    sign_tt(k)=sss/(xfpi(2));
    tau_tt(k)=ttt/(xfpi(2));
end

sign_tt;
tau_tt;

plot_graph=0;
if plot_graph==1
    loglog(dfs,abs(-sign_tt),'-.r')
    hold on
    loglog(dfs,abs(-sign_ttt),'-.k')
    %xlabel('distance sill-faille (m)')
    %ylabel('taux de variation de \sigma_n par Pa de surpression dans le sill')
    %title('variation de la contrainte normale en fonction de la distance sill-faille','fontsize',15)

    loglog(dfs,abs(tau_tt),'k')
    hold on
    loglog(dfs,abs(tau_ttt),'r')
    legend('\sigma_n (centre)','\sigma_n (mean)', '\tau (centre)','\tau (mean)','fontsize',22)
    xlabel('distance sill-faille (m)','fontsize',22)
    if mag_lens_act==1
        ylabel('taux de variation des differentes contraintes par Pa de surpression dans le sill','fontsize',18)
    elseif dyke_act==1
        ylabel("taux de variation des differentes contraintes par m d'ouverture du dike",'fontsize',22)
    else 
        ylabel('??','fontsize',22)
    end
    title('variation absolue des contraintes en fonction de la distance sill-faille','fontsize',25)
end

ploti=0;
if ploti==1
    h=pcolor(x,z,(SIGT.xx(:,:)));
    %h=pcolor(x,z,((1/2)*(SIGT.xx(:,:)+SIGT.zz(:,:)+2*(SIGT.xz(:,:)))));
    %h=pcolor(x,y,(sig_xx(:,:,1)));
    %h=pcolor(x,y,(sig_xy(:,:,1)));
    %h=pcolor(x,z,(sign_tot(:,:)));
    %h=pcolor(x,y,(tau_tot(:,:,2)-tau_tot(:,:,1))/dp);
    %h=pcolor(x,z,(tau_tot(:,:)));
    %h=pcolor(x,z,(SIGT.xx(:,:)+SIGT.zz(:,:)+SIGT.xz(:,:)));
    %h=pcolor(x,z,SIG_faultbot.xx(:,:));
    %h=pcolor(x,z,SIGT.zz(:,:));
    axis equal
    colormap(redblue)
    caxis([-1e-1,1e-1])
    colorbar
    hold on
    p=plot(xfp,zfp,'k');
    p.LineWidth=2;
    set(h,'EdgeColor','none')
    set(gca,'YDir')%,'reverse')
    ylabel('Profondeur (m)','fontsize',18)
    xlabel('distance au milieu du sill (m)','fontsize',18)
    title('$\sig_{xx}$','fontsize',22)
end

%plot(x,U2(1101,:))

dsigmwd=-sign_ttt_wd;
%dsigmwd.t=tau_ttt_wd;
dsigm.n=-sign_tt;
dsigm.t=tau_tt;

dsigOP=[];
for l=1:sxop(2)
    dsigOP=[dsigOP,SIGT.xx(izOP(1),ixOP(l))];
end
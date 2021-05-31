%--------------%
%--parameters--%
%--------------%
code_patin_ressort=0;

year=3600*24*365;  %--number of seconds in one year

if code_patin_ressort==1
   'interaction activated'
else
    mag_lens_act=1; 
    lithos_act=0;
    dyke_act=1;
    fault_act=0;
end

%--INTERACTIONS PARAMETERS--%
dfs=1000;%[0,1000,2000,3000,4000,5000,10000,15000,20000]; % distance btw end of the fault and sill (m)

%--SILL--%
if code_patin_ressort==1
   'interaction activated';
else
   P=0;    % overpressure initial in the crack (Pa)
   dp=1e8/year; % overpressure rate (Pa/s)
end
L=502;    % crack width (m) ; take an even nbr
theta=0; % crack angle relative to the x axis (deg) non implementé
y0=-1000;   % crack depth (m)
x0=0; % crack horizontal location (m)
nus=0.254;   % Poisson coeff

%--LITHOSTATIC PRESSURE--%
g=9.8; %--free-fall cte (m/s)
rho=2900; %--rock density (kg/m3)

%--DYKE--%
dd=0; % depth to top of crack (m)
Ld=y0; % height of crack (m)
b=0; % crack opening (m)
db=1;
mu=10e9; % shear modulus (Pa)
nud=nus; % poisson coefficient
A=1;%%mu*b/(2*pi*(1-nud)); 

sign_tt=dfs;
tau_tt=dfs;
sign_ttt=dfs;
tau_ttt=dfs;
rr=size(dfs);
SIG_lithos={};
SIG_dyke={};
SIG_fault={};
SIG_Mag={};
SIGT={};

%--BOX--%
pas=20;

%--FAULT--%
thetaf=45;%+90;%+180? % fault dip
Lf=abs(y0)/sind(thetaf); % surface of displacement on fault (?)
nu=nud;


%Stress-pertubation above a magma lense approximated

for k=1:rr(2)
    %--BOX--%
    d_fs=dfs(k);
    box=[dfs(k)+15000,abs(y0)+10000];
    x =-box(1):pas:box(1);
    y = -box(2):pas:0;
    tps=1:2;
    [XX,YY,TT] = meshgrid(x,y,tps);
    [X,Y] = meshgrid(x,y);
    s=size(XX);
    tau_tot=XX;
    sign_tot=XX;
   
    xsi=[x0+L/2+d_fs+(Lf*cosd(thetaf))/2,y0+(Lf*sind(thetaf))/2]; % origin dislocation
    df=xsi(2)+(Lf*sind(thetaf))/2; % fault top depth 
    if cosd(thetaf)==0
        xfp=zeros([1 10]);
        yfp=xsi(2)-(Lf*sind(thetaf))/2:(Lf*sind(thetaf))/9:xsi(2)+(Lf*sind(thetaf))/2;
    elseif sind(thetaf)==0
        xfp=xsi(1)-(Lf*cosd(thetaf))/2:(Lf*cosd(thetaf))/9:xsi(1)+(Lf*cosd(thetaf))/2;
        yfp=zeros([1 10]);
    else
        xfp=xsi(1)-(Lf*cosd(thetaf))/2:(Lf*cosd(thetaf))/10:xsi(1)+(Lf*cosd(thetaf))/2;
        yfp=xsi(2)-(Lf*sind(thetaf))/2:(Lf*sind(thetaf))/10:xsi(2)+(Lf*sind(thetaf))/2;
    end
    % where is le milieu de la faille ?
    F=X>xsi(1);
    Arrr=X(F);
    idxx=find(X(1,:,1)==Arrr(1));
    F=Y>xsi(2);
    Arrr=Y(F);
    idxy=find(Y(:,1,1)==Arrr(1));
    % where is la moyenne de la faille ?
    xfpi=size(xfp);
    idxfp=xfp;
    idyfp=yfp;
    for l=1:xfpi(2)
        idxfp(l)=find(X(1,:,1)>xfp(l)-pas/2 & X(1,:,1)<xfp(l)+pas/2);
        idyfp(l)=find(Y(:,1,1)>yfp(l)-pas/2 & Y(:,1,1)<yfp(l)+pas/2);
    end
    for tps=1:s(3)
        if mag_lens_act==1
            P=P+dp;
        elseif dyke_act==1
            b=b+db;
        else
            abc=2;
        end
        if lithos_act==1   
            SIG_lithos=lithos(rho,g,y);
        else
            SIG_lithos.xx = zeros(size(X));
            SIG_lithos.yy = zeros(size(X));
            SIG_lithos.xy = zeros(size(X));
        end
        if mag_lens_act==1
            SIG_Mag = mag_lens_f(X,Y,P,L,nus,y0,x0);
        else
            SIG_Mag.xx = zeros(size(X));
            SIG_Mag.yy = zeros(size(X));
            SIG_Mag.xy = zeros(size(X));
        end
        if dyke_act==1
            SIG_dyke = dyke(dd,Ld,X,Y,b,nud,mu);
        else
            SIG_dyke.xx = zeros(size(X));
            SIG_dyke.yy = zeros(size(X));
            SIG_dyke.xy = zeros(size(X));
        end
        if fault_act==1
            [SIG_fault,U1,U2] = fault_okada(X,Y,mu,Lf,df,nus,thetaf,xsi);
        else
            SIG_fault.xx = zeros(size(X));
            SIG_fault.yy = zeros(size(X));
            SIG_fault.xy = zeros(size(X));
            U1=0;
            U2=0;
        end
        SIGT.xx = SIG_Mag.xx+SIG_lithos.xx/abs(-rho*g)+SIG_dyke.xx/abs(A)+SIG_fault.xx/(mu*abs(Lf)/(2*pi*(1-nus)));
        SIGT.yy = SIG_Mag.yy+SIG_lithos.yy/abs(-rho*g)+SIG_dyke.yy/abs(A)+SIG_fault.yy/(mu*abs(Lf)/(2*pi*(1-nus)));
        SIGT.xy = SIG_Mag.xy+SIG_lithos.xy/abs(-rho*g)+SIG_dyke.xy/abs(A)+SIG_fault.xy/(mu*abs(Lf)/(2*pi*(1-nus)));
        tau_tot(:,:,tps)=(SIGT(1).xx-SIGT(1).yy)*cosd(thetaf)*sind(thetaf)+SIGT(1).xy*(sind(thetaf)^2-cosd(thetaf)^2); % repère de pierre
        sign_tot(:,:,tps)=SIGT(1).xx*sind(thetaf)^2+SIGT(1).yy*cosd(thetaf)^2-2*SIGT(1).xy*cosd(thetaf)*sind(thetaf); % repère de pierre
      
    end   
    sign_ttt(k)=(sign_tot(idxy,idxx,2)-sign_tot(idxy,idxx,1))/tps;
    tau_ttt(k)=(tau_tot(idxy,idxx,2)-tau_tot(idxy,idxx,1))/tps;
    sss=0;
    ttt=0;
    for m=2:xfpi(2)
        sss=sss+sign_tot(idyfp(m),idxfp(m),2)-sign_tot(idyfp(m),idxfp(m),1);
        ttt=ttt+tau_tot(idyfp(m),idxfp(m),2)-tau_tot(idyfp(m),idxfp(m),1);
    end
    sign_tt(k)=sss/(tps*xfpi(2));
    tau_tt(k)=ttt/(tps*xfpi(2));
    if mag_lens_act==1
        P=0;
    elseif dyke_act==1
        b=0;
    else
        abc=2;
    end
end

sign_tt
tau_tt

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
    h=pcolor(x,y,(SIGT.yy(:,:)+SIGT.xx(:,:)));
    %h=pcolor(x,y,(sig_xx(:,:,1)));
    %h=pcolor(x,y,(sig_xy(:,:,1)));
    %h=pcolor(x,y,(sign_tot(:,:,2)));
    %h=pcolor(x,y,(tau_tot(:,:,2)-tau_tot(:,:,1))/dp);
    %h=pcolor(x,y,(tau_tot(:,:,1)));
    %h=pcolor(x,y,((1/2)*(SIGT.xx(:,:)+SIGT.yy(:,:)+2*(SIGT.xy(:,:)))));
    %h=pcolor(x,y,((1/2)*(SIGT.xx(:,:)+SIGT.yy(:,:)-2*(SIGT.xy(:,:)))));
    %h=pcolor(x,y,SIGT(1).xx(:,:));
    axis equal
    colormap(redblue)
    %caxis([-0.1e7,0.1e7])
    colorbar
    hold on
    p=plot(xfp,yfp,'k');
    p.LineWidth=2;
    set(h,'EdgeColor','none')
    set(gca,'YDir')%,'reverse')
    ylabel('Profondeur (m)','fontsize',18)
    xlabel('distance au milieu du sill (m)','fontsize',18)
    title('Trace du tenseur des contraintes','fontsize',22)
end

dsigm.n=sign_tt;
dsigm.t=tau_tt;


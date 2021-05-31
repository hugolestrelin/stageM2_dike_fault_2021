%Stress-pertubation above a magma lense approximated
% tester avec theta recalculer à chaque cadrant
%--------------%
%--parameters--%
%--------------%

mag_lens_act=1;
lithos_act=0;
dyke_act=0;
fault_act=0;

%--INTERACTIONS PARAMETERS--%
dfs=[0,1000,2000,3000,4000,5000,10000,15000,20000]; % distance btw end of the fault and sill (m)

%--SILL--%
P=0;    % overpressure initial in the crack (Pa)
dp=-1; % overpressure rate (Pa/s)
L=2000;    % crack width (m)
theta=0; % crack angle relative to the x axis (deg) non implementé
y0=-3000;   % crack depth (m)
x0=0; % crack location (m)
nus=0.254;   % Poisson coeff

%--LITHOSTATIC PRESSURE--%
g=9.8; %--free-fall cte (m/s)
rho=2900; %--rock density (kg/m3)

%--DYKE--%
dd=0; % depth to top of crack (m)
Ld=-3000; % height of crack (m)
b=1; % crack opening (m)
mu=10e9; % shear modulus (Pa)
nud=nus; % poisson coefficient
A=1;%%mu*b/(2*pi*(1-nud)); 

%sig_ext=X;
%sig_comp=X;

%sig_zz=X;

%alpha_tot=X;
%U=X;
%V=X;
%U1_tot=X;
%U2_tot=X;

%SIG_tot = {};
%SIG_Mag = {};
%SIG_lithos = {};
%tau={};
%alpha={};

sign_tt=dfs;
tau_tt=dfs;
rr=size(dfs);

for k=1:rr(2)
    %--BOX--%
    d_fs=dfs(k);
    pas=200;
    box=[dfs(k)+5000,abs(y0)+1000];
    x =-box(1):pas:box(1);
    y = -box(2):pas:0;
    t=1:2;
    [X,Y,T] = meshgrid(x,y,t);
    s=size(X);
    
    tau_tot=X;
    sign_tot=X;
    sig_xx=X;
    sig_yy=X;
    sig_xy=X;
    %--FAULT--%
    thetaf=60;%+180? % fault dip
    Lf=abs(y0)/sind(thetaf); % surface of displacement on fault (?)
    xsi=[x0+L/2+d_fs+(Lf*cosd(thetaf))/2,y0+(Lf*sind(thetaf))/2]; % origin dislocation
    df=xsi(2)+(Lf*sind(thetaf))/2; % fault top depth 
    nu=nud;
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
    % where is la faille ?
    F=X>xsi(1);
    Arrr=X(F);
    idxx=find(X(1,:,1)==Arrr(1));
    F=Y>xsi(2);
    Arrr=Y(F);
    idxy=find(Y(:,1,1)==Arrr(1));
    for t=1:s(3)
        P=P+dp;
        if lithos_act==1   
            SIG_lithos=lithos(rho,g,y(i));
        else
            SIG_lithos=0;
        end
        if mag_lens_act==1
            SIG_Mag = mag_lens_f(X,Y,P,L,nus,y0,x0);
        else
            SIG_Mag = 0;
        end
        if dyke_act==1
            SIG_dyke = dyke(dd,Ld,X,Y,b,nud,mu);
        else
            SIG_dyke = 0;
        end
        if fault_act==1
            [SIG_fault,U1,U2] = fault_okada(X,Y,mu,Lf,df,nus,thetaf,xsi);
        else
            SIG_fault = 0;
            U1=0;
            U2=0;
        end
        SIGT = SIG_Mag+SIG_lithos/abs(-rho*g)+SIG_dyke/abs(A)+SIG_fault/(mu*abs(Lf)/(2*pi*(1-nus)));
        tau_tot(i,j,t)=(SIGT.xx-SIGT.yy)*cosd(thetaf)*sind(thetaf); 
        sign_tot(i,j,t)=SIGT.xx*cosd(thetaf)^2+SIGT.yy*sind(thetaf)^2-2*SIGT.xy*cosd(thetaf)*sind(thetaf);
    end   
    sign_tt(k)=(sign_tot(idxy,idxx,2)-sign_tot(idxy,idxx,1))/t;
    tau_tt(k)=(tau_tot(idxy,idxx,2)-tau_tot(idxy,idxx,1))/t;  
    P=0;
end

sign_tt
tau_tt

plot_graph=0;
if plot_graph==1
    plot(dfs,sign_tt)
    xlabel('distance sill-faille (m)')
    ylabel('sigma_n (Pa)')
    title('variation de la contrainte normale en fonction de la distance sill-faille','fontsize',15)

    plot(dfs,tau_tt)
    xlabel('distance sill-faille (m)')
    ylabel('tau (Pa)')
    title('variation de la contrainte tangentielle en fonction de la distance sill-faille','fontsize',15)
end

%displacement mag
%for i=1:s(1)
%    for j=1:s(2)
%        U(i,j)=U(i,j)*((sig_xx(i,j)+sig_yy(i,j))/2)/mean_normal_stress;
%        V(i,j)=V(i,j)*((sig_xx(i,j)+sig_yy(i,j))/2)/mean_normal_stress;
%    end
%end

ploti=1;
if ploti==1
    %h=pcolor(x,y,(sig_yy(:,:,1)+sig_xx(:,:,1))/2);
    %h=pcolor(x,y,(sig_xx(:,:,1)));
    %h=pcolor(x,y,(sig_xy(:,:,1)));
    h=pcolor(x,y,(sign_tot(:,:,2)-sign_tot(:,:,1))/dp);
    %h=pcolor(x,y,(tau_tot(:,:,2)-tau_tot(:,:,1))/dp);
    axis equal
    caxis([-1,1])
    colorbar
    hold on
    p=plot(xfp,yfp,'k');
    p.LineWidth=2;
    set(h,'EdgeColor','none')
    set(gca,'YDir')%,'reverse')
    title('mean stress','fontsize',15)
end

function [sig] = mag_lens_f(rx,ry,P,L,nu,wd,x0) % rx and ry matrix n*m*o ,tau,alpha]
    r=((rx-x0).^2+(ry-wd).^2).^(1/2);
    r1=((rx-x0+L/2).^2+(ry-wd).^2).^(1/2); %gauche <0
    r2=((L/2-rx-x0).^2+(ry-wd).^2).^(1/2); %droite >0
    theta=zeros(size(rx));
    theta1=zeros(size(rx));
    theta2=zeros(size(rx));
    theta(ry>wd & rx-x0==0)=rad2deg(pi/2);
    theta1(ry>wd & rx-x0==0)=atand(abs(ry(ry>wd & rx-x0==0)-wd)./abs(rx(ry>wd & rx-x0==0)-x0+L/2));
    theta2(ry>wd & rx-x0==0)=rad2deg(pi)-(atand(abs(ry(ry>wd & rx-x0==0)-wd)./abs(L/2-rx(ry>wd & rx-x0==0)-x0)));
    theta(ry>wd & (0<rx-x0)&&(rx-x0<L/2))=(atand(abs(ry(ry>wd & (0<rx-x0)&&(rx-x0<L/2))-wd)./abs(rx(ry>wd & (0<rx-x0)&&(rx-x0<L/2))-x0)));
    theta1(ry>wd & (0<rx-x0)&&(rx-x0<L/2))=(atand(abs(ry(ry>wd & (0<rx-x0)&&(rx-x0<L/2))-wd)./abs(rx(ry>wd & (0<rx-x0)&&(rx-x0<L/2))-x0+L/2)));
    theta2(ry>wd & (0<rx-x0)&&(rx-x0<L/2))=rad2deg(pi)-(atand(abs(ry(ry>wd & (0<rx-x0)&&(rx-x0<L/2))-wd)./abs(L/2-rx(ry>wd & (0<rx-x0)&&(rx-x0<L/2))-x0)));
    theta(ry>wd & rx-x0==L/2)=(atand(abs(ry(ry>wd & rx-x0==L/2)-wd)./abs(rx(ry>wd & rx-x0==L/2)-x0)));
    theta1(ry>wd & rx-x0==L/2)=(atand(abs(ry(ry>wd & rx-x0==L/2)-wd)./abs(rx(ry>wd & rx-x0==L/2)-x0+L/2)));
    theta2(ry>wd & rx-x0==L/2)=rad2deg((pi/2));
    theta(ry>wd & rx-x0>L/2)=(atand(abs(ry(ry>wd & rx-x0>L/2)-wd)./abs(rx(ry>wd & rx-x0>L/2)-x0)));
    theta1(ry>wd & rx-x0>L/2)=(atand(abs(ry(ry>wd & rx-x0>L/2)-wd)./abs(rx(ry>wd & rx-x0>L/2)-x0+L/2)));
    theta2(ry>wd & rx-x0>L/2)=(atand(abs(ry(ry>wd & rx-x0>L/2)-wd)./abs(L/2-rx(ry>wd & rx-x0>L/2)-x0)));
    theta(ry>wd & (-L/2<rx-x0)&&(rx-x0<0))=rad2deg(pi)-(atand(abs(ry(ry>wd & (-L/2<rx-x0)&&(rx-x0<0))-wd)./abs(rx(ry>wd & (-L/2<rx-x0)&&(rx-x0<0))-x0)));
    theta1(ry>wd & (-L/2<rx-x0)&&(rx-x0<0))=(atand(abs(ry(ry>wd & (-L/2<rx-x0)&&(rx-x0<0))-wd)./abs(rx(ry>wd & (-L/2<rx-x0)&&(rx-x0<0))+L/2-x0)));
    theta2(ry>wd & (-L/2<rx-x0)&&(rx-x0<0))=rad2deg(pi)-(atand(abs(ry(ry>wd & (-L/2<rx-x0)&&(rx-x0<0))-wd)./abs(L/2-(rx(ry>wd & (-L/2<rx-x0)&&(rx-x0<0))-x0))));
    theta(ry>wd & rx-x0==-L/2)=rad2deg(pi)-(atand(abs(ry(ry>wd & rx-x0==-L/2)-wd)./abs((rx(ry>wd & rx-x0==-L/2)-x0))));
    theta1(ry>wd & rx-x0==-L/2)=rad2deg(pi/2);
    theta2(ry>wd & rx-x0==-L/2)=rad2deg(pi)-(atand(abs(ry(ry>wd & rx-x0==-L/2)-wd)./abs(L/2-(rx(ry>wd & rx-x0==-L/2)-x0))));
    theta(ry>wd & rx-x0<-L/2)=rad2deg(pi)-(atand(abs(ry(ry>wd & rx-x0<-L/2)-wd)./abs(rx(ry>wd & rx-x0<-L/2)-x0)));
    theta1(ry>wd & rx-x0<-L/2)=rad2deg(pi)-(atand(abs(ry(ry>wd & rx-x0<-L/2)-wd)./abs(rx(ry>wd & rx-x0<-L/2)-x0+L/2)));
    theta2(ry>wd & rx-x0<-L/2)=rad2deg(pi)-(atand(abs(ry(ry>wd & rx-x0<-L/2)-wd)./abs(L/2-(rx(ry>wd & rx-x0<-L/2)-x0))));
    theta(ry==wd & rx-x0==0)=rad2deg(pi/2);
    theta2(ry==wd & rx-x0==0)=rad2deg(pi);
    theta2(ry==wd &(0<rx-x0)&&(rx-x0<L/2))=rad2deg(pi);
    theta2(ry==wd & rx-x0==L/2)=rad2deg(pi/2);
    theta(ry==wd & (-L/2<rx-x0)&&(rx-x0<0))=pi;
    theta2(ry==wd & (-L/2<rx-x0)&&(rx-x0<0))=rad2deg(pi);
    theta(ry==wd & rx-x0==-L/2)=rad2deg(pi);
    theta1(ry==wd & rx-x0==-L/2)=rad2deg(pi/2);
    theta2(ry==wd & rx-x0==-L/2)=rad2deg(pi);
    theta(ry==wd & rx-x0<-L/2)=rad2deg(pi);
    theta1(ry==wd & rx-x0<-L/2)=rad2deg(pi);
    theta2(ry==wd & rx-x0<-L/2)=rad2deg(pi);
    theta(ry<wd & rx-x0==0)=rad2deg(3*pi/2);
    theta1(ry<wd & rx-x0==0)=rad2deg(2*pi)-atand(abs(wd-ry(ry<wd & rx-x0==0))./abs(rx(ry<wd & rx-x0==0)-x0+L/2));
    theta2(ry<wd & rx-x0==0)=rad2deg(pi)+(atand(abs(wd-ry(ry<wd & rx-x0==0))./abs(L/2-rx(ry<wd & rx-x0==0)-x0)));
    theta(ry<wd & (0<rx-x0)&&(rx-x0<L/2))=rad2deg(2*pi)-(atand(abs(wd-ry(ry==wd & (0<rx-x0)&&(rx-x0<L/2)))./abs(rx(ry==wd & (0<rx-x0)&&(rx-x0<L/2))-x0)));
    theta1(ry<wd & (0<rx-x0)&&(rx-x0<L/2))=rad2deg(2*pi)-(atand(abs(wd-ry(ry==wd & (0<rx-x0)&&(rx-x0<L/2)))./abs(rx(ry==wd & (0<rx-x0)&&(rx-x0<L/2))-x0+L/2)));
    theta2(ry<wd & (0<rx-x0)&&(rx-x0<L/2))=rad2deg(pi)+(atand(abs(wd-ry(ry==wd & (0<rx-x0)&&(rx-x0<L/2)))./abs(L/2-rx(ry==wd & (0<rx-x0)&&(rx-x0<L/2))-x0)));
    theta(ry<wd & rx-x0==L/2)=rad2deg(2*pi)-(atand(abs(wd-ry(ry<wd & rx-x0==L/2))./abs(rx(ry<wd & rx-x0==L/2)-x0)));
    theta1(ry<wd & rx-x0==L/2)=rad2deg(2*pi)-(atand(abs(wd-ry(ry<wd & rx-x0==L/2))./abs(rx(ry<wd & rx-x0==L/2)-x0+L/2)));
    theta2(ry<wd & rx-x0==L/2)=rad2deg(3*pi/2);
    theta(ry<wd & rx-x0>L/2)=rad2deg(2*pi)-(atand(abs(wd-ry(ry<wd & rx-x0>L/2))./abs(rx(ry<wd & rx-x0>L/2)-x0)));
    theta1(ry<wd & rx-x0>L/2)=rad2deg(2*pi)-(atand(abs(wd-ry(ry<wd & rx-x0>L/2))./abs(rx(ry<wd & rx-x0>L/2)-x0+L/2)));
    theta2(ry<wd & rx-x0>L/2)=rad2deg(2*pi)-(atand(abs(wd-ry(ry<wd & rx-x0>L/2))./abs(L/2-rx(ry<wd & rx-x0>L/2)-x0)));
    theta(ry<wd & (-L/2<rx-x0)&&(rx-x0<0))=rad2deg(pi)+(atand(abs(wd-ry(ry<wd & (-L/2<rx-x0)&&(rx-x0<0)))./abs(rx(ry<wd & (-L/2<rx-x0)&&(rx-x0<0))-x0)));
    theta1(ry<wd & (-L/2<rx-x0)&&(rx-x0<0))=rad2deg(2*pi)-(atand(abs(wd-ry(ry<wd & (-L/2<rx-x0)&&(rx-x0<0)))./abs(rx(ry<wd & (-L/2<rx-x0)&&(rx-x0<0))-x0+L/2)));
    theta2(ry<wd & (-L/2<rx-x0)&&(rx-x0<0))=rad2deg(pi)+(atand(abs(wd-ry(ry<wd & (-L/2<rx-x0)&&(rx-x0<0)))./abs(L/2-(rx(ry<wd & (-L/2<rx-x0)&&(rx-x0<0))-x0))));
    theta(ry<wd & rx-x0==-L/2)=rad2deg(pi)+(atand(abs(wd-ry(ry<wd & rx-x0==-L/2))./abs(rx(ry<wd & rx-x0==-L/2)-x0)));
    theta1(ry<wd & rx-x0==-L/2)=rad2deg(3*pi/2);
    theta2(ry<wd & rx-x0==-L/2)=rad2deg(pi)+(atand(abs(wd-ry(ry<wd & rx-x0==-L/2))./abs(L/2-(rx(ry<wd & rx-x0==-L/2)-x0))));
    theta(ry<wd & rx-x0<-L/2)=rad2deg(pi)+(atand(abs(wd-ry(ry<wd & rx-x0<-L/2))./abs(rx(ry<wd & rx-x0<-L/2)-x0)));
    theta1(ry<wd & rx-x0<-L/2)=rad2deg(pi)+(atand(abs(wd-ry(ry<wd & rx-x0<-L/2))./abs(rx(ry<wd & rx-x0<-L/2)-x0+L/2)));
    theta2(ry<wd & rx-x0<-L/2)=rad2deg(pi)+(atand(abs(wd-ry(ry<wd & rx-x0<-L/2))./abs(L/2-(rx(ry<wd & rx-x0<-L/2)-x0))));
%     if ry>wd
%         if rx-x0==0
%             theta=rad2deg(pi/2);
%             theta1=atand(abs(ry-wd)/abs(rx-x0+L/2));
%             theta2=rad2deg(pi)-(atand(abs(ry-wd)/abs(L/2-rx-x0)));
%         elseif (0<rx-x0)&&(rx-x0<L/2)
%             theta=(atand(abs(ry-wd)/abs(rx-x0)));
%             theta1=(atand(abs(ry-wd)/abs(rx-x0+L/2)));
%             theta2=rad2deg(pi)-(atand(abs(ry-wd)/abs(L/2-rx-x0)));
%         elseif rx-x0==L/2
%             theta=(atand(abs(ry-wd)/abs(rx-x0)));
%             theta1=(atand(abs(ry-wd)/abs(rx-x0+L/2)));
%             theta2=rad2deg((pi/2));
%         elseif rx-x0>L/2
%             theta=(atand(abs(ry-wd)/abs(rx-x0)));
%             theta1=(atand(abs(ry-wd)/abs(rx-x0+L/2)));
%             theta2=(atand(abs(ry-wd)/abs(L/2-rx-x0)));
%         elseif (-L/2<rx-x0)&&(rx-x0<0)
%             theta=rad2deg(pi)-(atand(abs(ry-wd)/abs(rx-x0)));
%             theta1=(atand(abs(ry-wd)/abs(rx+L/2-x0)));
%             theta2=rad2deg(pi)-(atand(abs(ry-wd)/abs(L/2-(rx-x0))));
%         elseif rx-x0==-L/2
%             theta=rad2deg(pi)-(atand(abs(ry-wd)/abs((rx-x0))));
%             theta1=rad2deg(pi/2);
%             theta2=rad2deg(pi)-(atand(abs(ry-wd)/abs(L/2-(rx-x0))));
%         elseif rx-x0<-L/2
%             theta=rad2deg(pi)-(atand(abs(ry-wd)/abs(rx-x0)));
%             theta1=rad2deg(pi)-(atand(abs(ry-wd)/abs(rx-x0+L/2)));
%             theta2=rad2deg(pi)-(atand(abs(ry-wd)/abs(L/2-(rx-x0))));
%         end
%     elseif ry==wd
%         if rx-x0==0
%             theta=rad2deg(pi/2);
%             theta1=0;
%             theta2=rad2deg(pi);
%         elseif (0<rx-x0)&&(rx-x0<L/2)
%             theta=0;
%             theta1=0;
%             theta2=rad2deg(pi);
%         elseif rx-x0==L/2
%             theta=0;
%             theta1=0;
%             theta2=rad2deg(pi/2);
%         elseif rx-x0>L/2
%             theta=0;
%             theta1=0;
%             theta2=0;
%         elseif (-L/2<rx-x0)&&(rx-x0<0)
%             theta=pi;
%             theta1=0;
%             theta2=rad2deg(pi);
%         elseif rx-x0==-L/2
%             theta=rad2deg(pi);
%             theta1=rad2deg(pi/2);
%             theta2=rad2deg(pi);
%         elseif rx-x0<-L/2
%             theta=rad2deg(pi);
%             theta1=rad2deg(pi);
%             theta2=rad2deg(pi);
%         end
%     else
%         if rx-x0==0
%             theta=rad2deg(3*pi/2);
%             theta1=rad2deg(2*pi)-atand(abs(wd-ry)/abs(rx-x0+L/2));
%             theta2=rad2deg(pi)+(atand(abs(wd-ry)/abs(L/2-rx-x0)));
%         elseif (0<rx-x0)&&(rx-x0<L/2)
%             theta=rad2deg(2*pi)-(atand(abs(wd-ry)/abs(rx-x0)));
%             theta1=rad2deg(2*pi)-(atand(abs(wd-ry)/abs(rx-x0+L/2)));
%             theta2=rad2deg(pi)+(atand(abs(wd-ry)/abs(L/2-rx-x0)));
%         elseif rx-x0==L/2
%             theta=rad2deg(2*pi)-(atand(abs(wd-ry)/abs(rx-x0)));
%             theta1=rad2deg(2*pi)-(atand(abs(wd-ry)/abs(rx-x0+L/2)));
%             theta2=rad2deg(3*pi/2);
%         elseif rx-x0>L/2
%             theta=rad2deg(2*pi)-(atand(abs(wd-ry)/abs(rx-x0)));
%             theta1=rad2deg(2*pi)-(atand(abs(wd-ry)/abs(rx-x0+L/2)));
%             theta2=rad2deg(2*pi)-(atand(abs(wd-ry)/abs(L/2-rx-x0)));
%         elseif (-L/2<rx-x0)&&(rx-x0<0)
%             theta=rad2deg(pi)+(atand(abs(wd-ry)/abs(rx-x0)));
%             theta1=rad2deg(2*pi)-(atand(abs(wd-ry)/abs(rx-x0+L/2)));
%             theta2=rad2deg(pi)+(atand(abs(wd-ry)/abs(L/2-(rx-x0))));
%         elseif rx-x0==-L/2
%             theta=rad2deg(pi)+(atand(abs(wd-ry)/abs(rx-x0)));
%             theta1=rad2deg(3*pi/2);
%             theta2=rad2deg(pi)+(atand(abs(wd-ry)/abs(L/2-(rx-x0))));
%         elseif rx-x0<-L/2
%             theta=rad2deg(pi)+(atand(abs(wd-ry)/abs(rx-x0)));
%             theta1=rad2deg(pi)+(atand(abs(wd-ry)/abs(rx-x0+L/2)));
%             theta2=rad2deg(pi)+(atand(abs(wd-ry)/abs(L/2-(rx-x0))));
%         end
%     end
    a=L/2;
    sig.xx=P*(r./(r1.*r2).^(1/2).*(cosd(theta-(theta1+theta2)/2))-1-((a^2)*r./((r1.*r2).^(3/2))).*(sind(theta)).*(sind(3/2*(theta1+theta2))));
    sig.xy=P*(a^2)*r./(r1.*r2).^(3/2).*(sind(theta)).*(cosd(3/2*(theta1+theta2)));
    sig.yy=P*(r./(r1.*r2).^(1/2).*(cosd(theta-(theta1+theta2)/2))-1+((a^2)*r./((r1.*r2).^(3/2))).*(sind(theta)).*(sind(3/2*(theta1+theta2))));
    sig.zz=nu*(sigxx+sigyy);
    %tau=1/2*sqrt((sigxx-sigyy)^2*4*sigxy^2);
    %alpha=1/2*atand(2*sigxy/(sigxx-sigyy));
end

function [sig] = lithos(rho,g,z)

    sig.xx=-rho*g*z;
    sig.yy=-rho*g*z;
    sig.zz=-rho*g*z;
end

function [sig] = dyke(d,Ld,x,z,b,nu,mu) % normalisation sur box : d=d/y0 et Ld=Ld/L
%     % horizontal displacement
%     fxtop=disp_x_partial(d,x,z,b,nu);
%     fxbot=disp_x_partial(d+Ld,x,z,b,nu);
%     ux=fxbot-fxtop;
%     % vertical displacement
%     fztop=disp_z_partial(d,x,z,b,nu);
%     fzbot=disp_z_partial(d+Ld,x,z,b,nu);
%     uz=fzbot-fztop;
    % tau_xx
    ftxxtop=disp_txx_partial(d,x,z,b,nu,mu);
    ftxxbot=disp_txx_partial(d+Ld,x,z,b,nu,mu);
    txx=ftxxbot-ftxxtop;
    % tau_zz
    ftzztop=disp_tzz_partial(d,x,z,b,nu,mu);
    ftzzbot=disp_tzz_partial(d+Ld,x,z,b,nu,mu);
    tzz=ftzzbot-ftzztop;
    % tau_xz
    ftxztop=disp_txz_partial(d,x,z,b,nu,mu);
    ftxzbot=disp_txz_partial(d+Ld,x,z,b,nu,mu);
    txz=ftxzbot-ftxztop;
    sig.xx=txx;
    sig.xy=txz;
    sig.yy=tzz;
end

function [sig,u1,u2] = fault_okada(x,z,mu,L,d,nu,theta,xsi)
    s=[L*cosd(theta),L*sind(theta)];
    r1=((x-xsi(1)).^2+(z-xsi(2)).^2).^(1/2);
    r2=((x-xsi(1)).^2+(z+xsi(2)).^2).^(1/2);
    t1=atand((x-2*xsi(1))./(z-xsi(2)));
    t2=atand((x-2*xsi(1))./(z+xsi(2)));
    u1=s(1)/(pi*(1-nu))*((1-nu)/2*(t2-t1)+(x-xsi(1)).*(z-xsi(2))./(4*r1.^2)-(x-xsi(1)).*(z+(3-4*nu)*xsi(2))./(4*r2.^2)+xsi(2)*z.*(x-xsi(1)).*(z+xsi(2))./(r2.^4))+s(2)/(pi*(1-nu))*((1-2*nu)/4*log(r2./r1)-((z-xsi(2)).^2)./(4*r1.^2)+(z.^2+xsi(2)^2-4*(1-nu)*xsi(2)*(z+xsi(2)))./(4*r2.^2)+(xsi(2)*z.*(z+xsi(2)).^2)./(r2.^4));
    u2=s(1)/(pi*(1-nu))*((1-2*nu)/4*log(r2./r1)+((z-xsi(2)).^2)./(4*r1.^2)-((z+xsi(2)).^2-2*xsi(2)^2-2*(1-2*nu)*xsi(2)*(z+xsi(2)))./(4*r2.^2)+(xsi(2)*z.*(z+xsi(2)).^2)./(r2.^4))+s(2)/(pi*(1-nu))*((1-nu)/2*(t2-t1)+(x-xsi(1)).*(z-xsi(2))./(4*r1.^2)-(x-xsi(1)).*(z+(3-4*nu)*xsi(2))./(4*r2.^2)-xsi(2)*z.*(x-xsi(1)).*(z+xsi(2))./(r2.^4));
    sig.xx=mu*s(2)/(2*pi*(1-nu))*((x-xsi(1)).*((z-xsi(2)).^2-(x-xsi(1)).^2)./(r1.^4)-(x-xsi(1)).*((z+xsi(2)).^2-(x-xsi(1)).^2)./(r2.^4)+4*xsi(2)*(x-xsi(1))./(r2.^6)*((2*xsi(2)-z).*(z+xsi(2)).^2+(3*z+2*xsi(2)).*(x-xsi(1)).^2))-mu*s(1)/(2*pi*(1-nu))*((z-xsi(2)).*((z-xsi(2)).^2+3*(x-xsi(1)).^2)./(r1.^4)-(z+xsi(2)).*((z+xsi(2)).^2+3*(x-xsi(1)).^2)./(r2.^4)+2*xsi(2)./(r2.^6)*(6*z.*(z+xsi(2)).*(x-xsi(1)).^2-(z-xsi(2)).*(z+xsi(2)).^3-(x-xsi(1)).^4)); 
    sig.yy=-mu*s(2)/(2*pi*(1-nu))*(d*(x-xsi(1)).*(3*(z-xsi(2)).^2-(x-xsi(1)).^2)./(r1.^4)-(x-xsi(1)).*(3*(z+xsi(2)).^2+(x-xsi(1)).^2)./(r2.^4)-4*xsi(2)*(x-xsi(1)).*z./(r2.^6).*(3*(z+xsi(2)).^2-(x-xsi(1))))-mu*s(1)/(2*pi*(1-nu))*((z-xsi(2)).*((z-xsi(2)).^2-(x-xsi(1)).^2)./(r1.^4)-(z+xsi(2)).*((z+xsi(2)).^2-(x-xsi(1)).^2)./(r2.^4)-2*xsi(2)/(r2.^6).*(6*z*(z+xsi(2)).*(x-xsi(1)).^2-(3*z+xsi(2)).*(z+xsi(2)).^3+(x-xsi(1)).^4)); 
    sig.xy=mu*s(2)/(2*pi*(1-nu))*((z-xsi(2)).*((z-xsi(2)).^2-(x-xsi(1)).^2)./(r1.^4)-(z+xsi(2)).*((z+xsi(2)).^2-(x-xsi(1)).^2)./(r2.^4)+2*xsi(2)/(r2.^6)*(6*z.*(z+xsi(2)).*(x-xsi(1)).^2-(x-xsi(1)).^4+(xsi(2)-z).*(z+xsi(2)).^3))-mu*s(1)/(2*pi*(1-nu))*((x-xsi(1)).*((z-xsi(2)).^2-(x-xsi(1)).^2)./(r1.^4)-(x-xsi(1)).*((z+xsi(2)).^2-(x-xsi(1)).^2)./(r2.^4)+4*xsi(2)*z.*(x-xsi(1))./(r2.^6).*(3*(z+xsi(2)).^2-(x-xsi(1)).^2)); 
end
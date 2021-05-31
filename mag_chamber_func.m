% parameters 

a=1; % semimajor axis length
c=0.9*a; % distance btween center and focal points 
b=sqrt(a^2-c^2); % semiminor axis length % infini
nu=0.4; % poisson's ratio
mu=1; %module de cisaillement 
lam=mu; % coeff lamé
P=-0.1; 
chi=0; %???
y0=4000; %origine en y
z0=1; %origine en z
xsiv=[1,1,1];
theta=0; %deg

%espace
pas=100;
box=4000;
x =0;%-box:pas:box;
y =0:pas:2*box;
z = 0:pas:box;
[X,Y,Z] = meshgrid(x,y,z);
U1s=X;
U2s=X;
U3s=X;

U1={};
U2={};
U3={};
s=size(X);
for i=1:s(1)
    for j=1:s(2)
        for k=1:s(3)
            [U1,U2,U3]=mag_chamber(xsiv,X(i,j,k),Y(i,j,k),Z(i,j,k),y0,z0,a,b,c,nu,mu,P,chi,lam,theta);
            U1s(i,j,k)=U1;
            U2s(i,j,k)=U2;
            U3s(i,j,k)=U3;
        end
    end
end

[A,B] = meshgrid(y,z);
U1d=A;
U2d=A;
U3d=A;
for j=1:s(1)
    for k=1:s(3)
            U1d(k,j)=real(U1s(j,round(s(2)/2),k));
            U2d(k,j)=real(U2s(j,round(s(2)/2),k));
            U3d(k,j)=real(U3s(j,round(s(2)/2),k));
   end
end

%imagesc(U3d)
plot(A(1,:),-U3d(1,:));

function [U1,U2,U3] = mag_chamber(xsiv,x,y,z,y0,z0,a,b,c,nu,mu,P,chi,lam,theta)
    xsi=xsiv(1)+xsiv(2)+xsiv(3);
    Q=3/(8*pi*(1-nu));
    R=(1-2*nu)/(8*pi*(1-nu));
    Ia=-2*pi*a*b^2*(2/(a*c^2)+1/(c^3)*log((a-c)/(a+c)));
    Iaa=-2*pi*a*b^2*(2/(3*a^3*c^2)+2/(a*c^4)+1/c^5*log((a-c)/(a+c)));
    a11=Q*a^2*Iaa+R*Ia-1;
    a12=(Q-R)*Ia-1;
    a21=2*R*(Ia-4*pi);
    a22=-16*pi*R;
    B1=1/(3*lam+2*mu)*(3*lam*(Q-R)*Ia-(3*lam-2*mu));
    B2=6/(3*lam+2*mu)*(mu-8*pi*lam*R);
    Pbarre=P*(B1*a22-B2*a12)/(a11*a22-a12*a21);
    Pstar=P*(B2*a11-B1*a21)/(a11*a22-a12*a21);
    a1=-2*b^2*Pbarre;
    b1=3*(b^2)/(c^2)*Pbarre+2*(1-2*nu)*Pstar;
    Pd=pi*a*b^2/(2*mu*c^3)*(c^2-chi^2)*Pbarre;
    Pc=pi*a*b^2/(2*mu*c^3)*(a1+b1*(c^2-chi^2));
    
    x1=x;
    x2=y-y0;
    x3=z-z0;
    x3b=z+z0;
    y1=x1;
    y2=x2-xsiv(2);
    y3=x3-xsiv(3);
    y3b=x3+xsiv(3);
    r2=x2*sind(theta)-x3*cosd(theta);
    r3=x2*cosd(theta)+x3*sind(theta);
    q2=x2*sind(theta)+x3b*cosd(theta);
    q3=-x2*cosd(theta)+x3b*sind(0);
    r3b=r3-xsi;
    q3b=q3+xsi;
    R1=sqrt(y1^2+y2^2+y3^2);
    R2=sqrt(y1^2+y2^2+y3b^(2));
    C0=y0*cosd(theta)+z0*sind(theta);
    beta=(cosd(theta)*q2+(1+sind(theta))*(R2+q3b))/(cosd(theta)*y1);
    Astar1=a1/(R1*(R1+r3b))+b1*(log(R1+r3b)+(r3+xsi)/(R1+r3b));
    Astar1b=-a1/(R2*(R2+q3b))-b1*(log(R2+q3b)+(q3-xsi)/(R2+q3b));
    A1=xsi/R1+log(R1+r3b);
    A1b=xsi/R2-log(R2+q3b);
    A2=R1-r3*log(R1+r3b);
    A2b=R2-q3*log(R2+q3b);
    A3=(xsi*r3b)/R1+R1;
    A3b=(xsi*q3b)/R2-R2;
    Bstar=a1/R1+2*b1*A2+(3-4*nu)*(a1/R2+2*b1*A2b);
    Fstar1=2*z*(cosd(theta)*q2*(a1*(2*R2+q3b)/(R2^3*(R2+q3b)^2)-b1*(R2+2*xsi)/(R2*(R2+q3b)^2))+sind(theta)*(a1/(R2^3)-2*b1*(R2+xsi)/(R2*(R2+q3b))));
    Fstar2=2*z*(a1*y3b/(R2^3)-2*b1*(sind(theta)*A1b+cosd(theta)*q2*(R2+xsi)/(R2*(R2+q3b))));
    B=xsi*(xsi+C0)/R2-A2b-C0*log(R2*q3b);
    F1=-2*sind(theta)*z*(xsi*(xsi+C0)/(R2^3)+(R2+xsi+C0)/(R2*(R2+q3b))+4*(1-nu)*(R2+xsi)/(R2*(R2+q3b)));
    F2=-2*sind(theta)*z*(xsi*(xsi+C0)*q3b/(R2^3)+C0/R2+(5-4*nu)*A1b);
    f1=xsi*y1/(R2+y3b)+3/(cosd(theta)^2)*(y1*sind(theta)*log(R2+y3b)-y1*log(R2+q3b)+2*q2*atand(beta))+2*y1*log(R2+q3b)-4*x3b/cosd(theta)*atand(beta);
    f2=xsi*y2/(R2+y3b)+(3/cosd(theta)^2)*(q2*sind(theta)*log(R2+q3b)-q2*log(R2+y3b)+2*y1*sind(theta)*atand(beta)+cosd(theta)*(R2-xsiv(3)))-2*cosd(theta)*A2b+2/cosd(theta)*(x3b*log(R2+y3b)-q3*log(R2+q3b));
    f3=1/cosd(theta)*(q2*log(R2+q3b)-q2*sind(theta)*log(R2+q3b)+2*y1*atand(beta))+2*sind(theta)*A2b+q3*log(R2+y3b)-xsi;
    
    Ustar1=1/(16*mu*(1-nu))*a*b^2/(c^3)*(Astar1*y1+(3-4*nu)*Astar1b*y1+Fstar1*y1);
    Ubarre1=1/(8*mu*(1-nu))*a*b^2/(c^3)*Pbarre*((A1*y1+(3-4*nu)*A1b*y1+F1*y1)-4*(1-nu)*(1-2*nu)*f1);
    Ustar2=1/(16*mu*(1-nu))*a*b^2/(c^3)*(sind(theta)*(Astar1*r2+(3-4*nu)*Astar1b*q2+Fstar1*q2)+cosd(theta)*(Bstar-Fstar2)+2*sind(theta)*cosd(theta)*z*Astar1b);
    Ubarre2=1/(8*mu*(1-nu))*a*b^2/(c^3)*Pbarre*(sind(theta)*(A1*r2+(3-4*nu)*A1b*q2+F1*q2)-4*(1-nu)*(1-2*nu)*f2+4*(1-nu)*cosd(theta)*(A2+A2b)+cosd(theta)*(A3-(3-4*nu)*A3b-F2));
    Ustar3=1/(16*mu*(1-nu))*a*b^2/(c^3)*(-cosd(theta)*(Astar1*r2+(3-4*nu)*Astar1b*q2-Fstar1*q2)+sind(0)*(Bstar+Fstar2)+2*cosd(theta)^2*z*Astar1b);
    Ubarre3=1/(8*mu*(1-nu))*a*b^2/(c^3)*Pbarre*(cosd(theta)*(-A1*r2+(3-4*nu)*A1b*q2+F1*q2)+4*(1-nu)*(1-2*nu)*f3+4*(1-nu)*sind(theta)*(A2+A2b)+sind(theta)*(A3+(3-4*nu)*A3b+F2-2*(3-4*nu)*B));
    
    U1=Ustar1+Ubarre1;
    U2=Ustar2+Ubarre2;
    U3=Ustar3+Ubarre3;
    
    %eps11=U1/dx1;
end

% solution en far-field ; 90° ?
function [Ur,U3] = mag_chamber_ff(xsiv,x,y,z,y0,z0,a,b,c,nu,mu,P,chi,lam,theta)
    xsi=xsiv(1)+xsiv(2)+xsiv(3);
    Q=3/(8*pi*(1-nu));
    R=(1-2*nu)/(8*pi*(1-nu));
    Ia=-2*pi*a*b^2*(2/(a*c^2)+1/(c^3)*log((a-c)/(a+c)));
    Iaa=-2*pi*a*b^2*(2/(3*a^3*c^2)+2/(a*c^4)+1/c^5*log((a-c)/(a+c)));
    a11=Q*a^2*Iaa+R*Ia-1;
    a12=(Q-R)*Ia-1;
    a21=2*R*(Ia-4*pi);
    a22=-16*pi*R;
    B1=1/(3*lam+2*mu)*(3*lam*(Q-R)*Ia-(3*lam-2*mu));
    B2=6/(3*lam+2*mu)*(mu-8*pi*lam*R);
    Pbarre=P*(B1*a22-B2*a12)/(a11*a22-a12*a21);
    Pstar=P*(B2*a11-B1*a21)/(a11*a22-a12*a21);
    a1=-2*b^2*Pbarre;
    b1=3*(b^2)/(c^2)*Pbarre+2*(1-2*nu)*Pstar;
    Pd=pi*a*b^2/(2*mu*c^3)*(c^2-chi^2)*Pbarre;
    Pc=pi*a*b^2/(2*mu*c^3)*(a1+b1*(c^2-chi^2));
    
    x1=x;
    x2=y-y0;
    x3=z-z0;
    x3b=z+z0;
    y1=x1;
    y2=x2-xsiv(2);
    y3=x3-xsiv(3);
    y3b=x3+xsiv(3);
    r2=x2*sind(theta)-x3*cosd(theta);
    r3=x2*cosd(theta)+x3*sind(theta);
    q2=x2*sind(theta)+x3b*cosd(theta);
    q3=-x2*cosd(theta)+x3b*sind(0);
    r3b=r3-xsi;
    q3b=q3+xsi;
    R1=sqrt(y1^2+y2^2+y3^2);
    R2=sqrt(y1^2+y2^2+y3^(-2));
    C0=y0*cosd(theta)+z0*sind(theta);
    beta=(cosd(theta)*q2+(1+sind(theta))*(R2+q3b))/(cosd(theta)*y1);
    
    Ar=r/(R1^3)-r/(R2^3);
    Az=x3/(R1^3)+x3b/(R2^3);
    Br=3*r*((x3^2)/(R1^5)-(x3b^2)/(R2^5));
    Bz=3*((x3^3)/(R1^5)+(x3b^3)/(R2^5));
    F1=4*(1-nu)*(2*nu/(R2^3)-3*x3b^2/(R2^5));
    F2=6*z*((z0+x3b)/(R2^5)-5*z0*(x3b^2)/(R2^7));
    Fc=4*(1-nu)/(R2^3)+6*z*x3b/(R2^5);
end
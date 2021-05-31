function f=disp_z_partial(s,x2,x3,b,nu)

del=pi/2;
R2=(x2-s*cos(del)).^2+(x3-s*sin(del)).^2;
S2=(x2-s*cos(del)).^2+(x3+s*sin(del)).^2;

A=b/(4*pi*(1-nu));
f=A*(2*(1-nu)*(log(sqrt(R2./S2))-2*s*(x3+s)./S2)-(log(sqrt(R2./S2))-x2.^2./R2-(x3.^2 +s.^2)./S2-4*s*x2.^2.*x3.*S2.^-2));
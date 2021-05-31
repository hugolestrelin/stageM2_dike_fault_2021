function f=disp_txz_partial(s,x2,x3,b,nu,mu)

del=pi/2;
R2=(x2-s*cos(del)).^2+(x3-s*sin(del)).^2;
S2=(x2-s*cos(del)).^2+(x3+s*sin(del)).^2;

A=mu*b/(2*pi*(1-nu));
f=-A.*x2.*(-(1./R2+1./S2)+2.*x2.^2.*R2.^-2+2*(x3.^2+s.^2-4*s.*x3).*S2.^-2+16*x2.^2.*x3.*s.*S2.^-3);
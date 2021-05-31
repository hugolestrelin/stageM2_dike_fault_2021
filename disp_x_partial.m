function f=disp_x_partial(s,x2,x3,b,nu)

del=pi/2;
R2=(x2-s*cos(del)).^2+(x3-s*sin(del)).^2;
S2=(x2-s*cos(del)).^2+(x3+s*sin(del)).^2;

A=b/(4*pi*(1-nu));

%f=A.*(2*(1-nu)*(2*x2.*s./S2+atan(x2./(x3-s))-atan(x2./(x3+s)))-(x2.*(x3-s).*(1./R2-1./S2)+4*x2.*x3.*s.*(x3+s).*S2.^-2));


a1=atan2(x2,(x3-s));
a2=atan2(x2,(x3+s));
f=A.*(2*(1-nu)*(2*x2.*s./S2+a1-a2)-(x2.*(x3-s).*(1./R2-1./S2)+4*x2.*x3.*s.*(x3+s).*S2.^-2));

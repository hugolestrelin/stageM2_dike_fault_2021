function [sig] = lithos(rho,g,z)
    sig.xx=-rho*g*z;
    sig.yy=-rho*g*z;
    sig.zz=-rho*g*z;
end
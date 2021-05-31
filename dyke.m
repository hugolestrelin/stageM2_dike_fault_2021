function [sig] = dyke(d,Ld,x,z,b,nu,mu) 
    ftxxtop=disp_txx_partial(-d,x,-z,b,nu,mu);
    ftxxbot=disp_txx_partial(-(d+Ld),x,-z,b,nu,mu);
    txx=ftxxbot-ftxxtop;
    ftzztop=disp_tzz_partial(-d,x,-z,b,nu,mu);
    ftzzbot=disp_tzz_partial(-(d+Ld),x,-z,b,nu,mu);
    tzz=ftzzbot-ftzztop;
    ftxztop=disp_txz_partial(-d,x,-z,b,nu,mu);
    ftxzbot=disp_txz_partial(-(d+Ld),x,-z,b,nu,mu);
    txz=ftxzbot-ftxztop;
    sig.xx=txx;
    sig.xz=txz;
    sig.zz=tzz;
end

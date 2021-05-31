clear all


% calculates the stress field and displacement associated 
% with a vertical mode-I dislocation in an elastic half-space.

% J.A. Olive, June 2012 after Singh and Singh, 2000


% COORDINATES (x,y) (choose even number of points to avoid singularities)
xx=linspace(-20e3,20e3,400);
zz=linspace(0,10e3,100);
[x z]=meshgrid(xx,zz);


% DISLOCATION see setup in Fig. 1, Singh and Singh 2000
d=1000; % depth to top of crack (m)
L=3000; % height of crack (m)
b=2; % crack opening (m)

% ELASTIC PARAMETERS
mu=30e9; % shear modulus (Pa)
nu=0.25; % poisson coefficient






% horizontal displacement

fxtop=disp_x_partial(d,x,z,b,nu);
fxbot=disp_x_partial(d+L,x,z,b,nu);
ux=fxbot-fxtop;

% vertical displacement

fztop=disp_z_partial(d,x,z,b,nu);
fzbot=disp_z_partial(d+L,x,z,b,nu);
uz=fzbot-fztop;

% tau_xx

ftxxtop=disp_txx_partial(d,x,z,b,nu,mu);
ftxxbot=disp_txx_partial(d+L,x,z,b,nu,mu);
txx=ftxxbot-ftxxtop;

% tau_zz

ftzztop=disp_tzz_partial(d,x,z,b,nu,mu);
ftzzbot=disp_tzz_partial(d+L,x,z,b,nu,mu);
tzz=ftzzbot-ftzztop;

% tau_xz

ftxztop=disp_txz_partial(d,x,z,b,nu,mu);
ftxzbot=disp_txz_partial(d+L,x,z,b,nu,mu);
txz=ftxzbot-ftxztop;



%PLOT
figure
nsmp=5;
h=pcolor(xx,zz,txx);
colorbar
set(h,'EdgeColor','none')
set(gca,'YDir','reverse')
axis equal
caxis([-1*10^8,1*10^8])
hold on
quiver(x(1:nsmp:end,1:nsmp:end),z(1:nsmp:end,1:nsmp:end),ux(1:nsmp:end,1:nsmp:end),uz(1:nsmp:end,1:nsmp:end))
title('\tau_{xx} (Pa)','fontsize',15)





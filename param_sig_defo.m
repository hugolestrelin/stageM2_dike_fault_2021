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
dp=1; % overpressure rate (Pa/s)
L=2000;    % crack width (m)
theta=0; % crack angle relative to the x axis (deg) non implementé
y0=-10000;   % crack depth (m)
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
thetaf=60;%+180? % fault dip
Lf=abs(y0)/sind(thetaf); % surface of displacement on fault (?)
nu=nud;
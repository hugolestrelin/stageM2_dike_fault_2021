%--------------%
%--parameters--%
%--------------%

mag_lens_act=0;
lithos_act=0;
dyke_act=1;
fault_act=0;

%--INTERACTIONS PARAMETERS--%
dfs=[1000];%[0,1000,2000,3000,4000,5000,10000,15000,20000]; % distance btw end of the fault and sill (m)

%--SILL--%
P=1;    % overpressure initial in the crack (Pa)
dp=1; % overpressure rate (Pa/s)
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
thetaf=45;%+180? % fault dip
Lf=abs(y0)/sind(thetaf); % surface of displacement on fault (?)
nu=nud;
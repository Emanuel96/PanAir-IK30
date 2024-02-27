% CONTROL POINTS
CPX = zeros([N,1]);
CPY = zeros([N,1]);
% VORTEX CORES
XM = zeros([NVs,1]);
YM = zeros([NVs,1]);
UM = zeros([NVs,1]);
VM = zeros([NVs,1]);
CV = zeros([NVs,1]);
% INFLUENCE MATRICES
%AN = zeros([N,N]);
%AT = zeros([N,N]);
A = zeros([N,N]);
BNI_W = zeros([N,1]);
BTI_W = zeros([N,1]);
AXW_J = zeros([N,1]);
AYW_J = zeros([N,1]);
BXW_J = zeros([N,1]);
BYW_J = zeros([N,1]);
AXH_J = zeros([NTS,N]);
AYH_J = zeros([NTS,N]);
BXH_J = zeros([NTS,N]);
BYH_J = zeros([NTS,N]);
BXH_W = zeros([NTS,1]);
BYH_W = zeros([NTS,1]);
CNI_M = zeros([N,NVs]);
CTI_M = zeros([N,NVs]);
CXW_M = zeros([1,NVs]);
CYW_M = zeros([1,NVs]);
CXH_M = zeros([NVs,NVs]);
CYH_M = zeros([NVs,NVs]);
% VELOCITIES VECTORS
PLUNG_VEL = zeros([NTS,1]);
ROT = zeros([NTS,1]);
ALPHA_1 = zeros([NTS,1]);
ALPHA_2 = zeros([NTS,1]);
ALPHA_3 = zeros([NTS,1]);
ALPHA_EFF = zeros([NTS,1]);

% NUMBER OF PANELS IN INTEGRATION PATH
N_PHI = 20;
% POSITION AND CONTROL POINTS OF PANELS
XZ = zeros([N_PHI+1,1]);
YZ = zeros([N_PHI+1,1]);
XMIDZ = zeros([N_PHI,1]);
YMIDZ = zeros([N_PHI,1]);
% VELOCITY AT PANELS
V_T_PHI_F = zeros([N_PHI,1]);
V_T_PHI = zeros([N,1]);
PHI_NODE = zeros([N+1,1]);
PHI = zeros([N,1]);
PHIM1 = zeros([N,1]);
% PRESSURE COEFFICIENT
%CP = zeros([N,1]);
% FORCES AND MOMENTS
CD = zeros([NTS,1]);
CL = zeros([NTS,1]);
CM = zeros([NTS,1]);
% TANGENTIAL VELOCITIES
V_t = zeros([N,1]);
% AIRFOIL VELOCITIES
V_XY_AIRFOIL = zeros([N,2]);
V_AL_AIRFOIL = zeros([N,2]);
V_AL_NODE_AIRFOIL = zeros([N+1,2]);
V_STREAM = zeros([N,2]);
% VECTORS FOR SOLVING SYSTEM OF EQUATIONS
b = zeros([N,1]);
c = zeros([N,1]);
% INITIALIZE TRAILING-EDGE PANEL PARAMETERS
UW_old = V_INF(1);
VW_old = V_INF(2);
UW = V_INF(1);
VW = V_INF(2);
Theta_old = atan2(VW,UW);
Delta_old = dt*sqrt(UW^2+VW^2);
Theta = atan2(VW,UW);
Delta = dt*sqrt(UW^2+VW^2);
% INITIALIZE AIRFOIL POSITION
HXO = 0;
HYO = 0;
% INITIALIZE FIRST WAKE PANEL
x(N+2) = x(N+1)+Delta*cos(Theta);
y(N+2) = y(N+1)+Delta*sin(Theta);

%% CALCULATIONS
% INITIALIZE ITERATION NUMBER
iter = 0;
% NUMBER OF CORE VORTICES
NM = 0;

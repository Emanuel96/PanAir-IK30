
%% SIMULATION PARAMETERS
vis = 1.7894e-5;
rho = 1.225;
% AERODYNAMIC CHORD
chord = 1;

% ANGLE OF ATTACK
meanAoA = [0]'*pi/180;
% REYNOLDS NUMBER
Re = [1e4]';

% NONDIMENSIONAL AMPLITUDE
NonDimAmpl = [0]';

% REDUCED FREQUENCY
%Strouhal = 0.1;
%RedFreq = pi*Strouhal./NonDimAmpl;

alpha_induced = 5;
RedFreq = tand(alpha_induced)./NonDimAmpl;
%RedFreq = [0.05 0.1 0.2 0.4 0.8 1.6]';

RedFreq = 1;

% CONDITION NUMBER
cond_no = 1;

% PITCHING AMPLITUDE
PitchAmp = -[5 10]'*pi/180;
% PHASE
phase = 90*pi/180;
% PIVOT (CHORD PERCENTAGE)
pivot = 0.3;
% ONLY LEADING EDGE?
onlyLE = 0; % 0 - STANDARD FLAPPING | 1 - LEADING-EDGE |
smoother = 117.7;
% ASYMMETRY OF KINEMATICS
zeta = 0.5;
% ANIMATION ?
animation = 1;
% SAVE DATA?
save_data = 0;

Z_0 = 1;
Z_1 = 1.1;

%AMP_AL_0 = -9*pi/180;

%old_CP = 0;

%% NUMERICAL PARAMETERS FOR PERIODIC PROBLEMS
% NUMBER OF SIMULATED PERIODS
NPs = 7;
% NUMBER OF TIME STEPS PER PERIOD
NTSPP = 100;
% NUMBER OF TIMESTEPS
NTS = NPs*NTSPP;
% MAX NUMBER OF ITERATIONS PER TIMESTEP
NITER = 100;
% MAX NUMBER OF VORTICES
NVs = 2*NTSPP;
% CONVERGENCE REPORT
reportConvergence = 5;

%% NUMERICAL PARAMETERS FOR OTHER PROBLEMS
% % NUMBER OF TIMESTEPS
% NTS = 1000;
% % TIMESTEP (seconds)
% dt = 0.1;
% % MAX NUMBER OF ITERATIONS PER TIMESTEP
% NITER = 100;
% % MAX NUMBER OF VORTICES
% NVs = 0.5*NTS;
% % CONVERGENCE REPORT
% reportConvergence = 5;

%% AIRFOIL DISCRETIZATION

% NUMBER OF BOUNDARY POINTS (ODD NUMBER)
N_P = 100;



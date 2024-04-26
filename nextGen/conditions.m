
%% SIMULATION PARAMETERS

% AIR
%vis = 1.7894e-5;
%rho = 1.225;
% WATER
vis = 0.001003;
rho = 998.2;

% AERODYNAMIC CHORD
chord = 0.20;
% ANGLE OF ATTACK
meanAoA = [0]'*pi/180;
% REYNOLDS NUMBER
Re = [1e5]';
% NONDIMENSIONAL AMPLITUDE
NonDimAmpl = [0.5]';
% REDUCED FREQUENCY
%Strouhal = 0.05;
%RedFreq = pi*Strouhal./NonDimAmpl;
alpha_induced = 10;
RedFreq = tand(alpha_induced)./NonDimAmpl;
%RedFreq = [0.05 0.1 0.2 0.4]';
% PITCHING AMPLITUDE
PitchAmp = -[5]'*pi/180;
%PitchAmp = -[0:5:alpha_induced]'*pi/180;
% PHASE
phase = 0*pi/180;
% PIVOT (CHORD PERCENTAGE)
pivot = 0.3;
% ONLY LEADING EDGE?
onlyLE = 0; % 0 - STANDARD FLAPPING | 1 - LEADING-EDGE |
smoother = 117.7;
% ASYMMETRY OF KINEMATICS
zeta = 0.5;

% LAMBDA (PITCHING NONSINUSOIDAL)
lambda = [0]';





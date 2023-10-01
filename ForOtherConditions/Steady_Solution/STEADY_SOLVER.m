%  UNSTEADY HSPM by Komasho (Based on TENG)
%  Last Update: MARCH 28th, 2022 by Emanuel Camacho

clc; close all; %clear all;

%% SIMULATION PARAMETERS

% AERODYNAMIC CHORD
chord = 1;
% ANGLE OF ATTACK
AoA = 0*pi/180;
% REYNOLDS NUMBER
Re = 1e6;
pivot = 0.3;

alpha0 = -0*pi/180;

for iRe = 1:size(Re)
    
    % AIRFOIL VELOCITY AND MAGNITUDE
    V_mag = Re(iRe)*1.7894e-5/(chord*1.225);
        
    %% NUMERICAL PARAMETERS
    
    %% AIRFOIL COORDINATES
    
    % ALERT: (COORDINATES MUST ALWAYS BE CREATED/READ ON THE COUNTERCLOCKWISE DIRECTION)
    % NUMBER OF BOUNDARY POINTS (ODD NUMBER)
    N_P = 200;
    if rem(N_P,2) == 0 % FORCE ODD NUMBER
        N_P = N_P + 1;
    end
    % NUMBER OF PANELS
    N = N_P-1;
    
    %% NACA0012 COORDINATE GENERATOR (WITH COSSINE SPACING)
    
    beta = linspace(0,2*pi,N_P)';
    x1 = 0.5*(1-cos(beta+pi));
    y1 = 0.594689181*(0.298222773*sqrt(x1)-0.127125232.*x1-0.357907906*x1.^2+0.291984971*x1.^3-0.105174606*x1.^4);
    x = (x1-pivot)*chord;
    y=([y1(1:(N_P-1)/2);-y1((N_P-1)/2+1:N_P)])*chord;
    y=flip(y);
    x=round(x,7);
    y=round(y,7);
    
    x_wd = x;
    y_wd = y;
    
    x_wd = x*cos(alpha0)-y*sin(alpha0);
    y_wd = x*sin(alpha0)+y*cos(alpha0);

    rotate = 0*pi/180;
    for i = 1:N
        eff_rot = 1/(1+exp(50*(x_wd(i))))*rotate;
        x(i) = x_wd(i)*cos(eff_rot)-y_wd(i)*sin(eff_rot);
        y(i) = x_wd(i)*sin(eff_rot)+y_wd(i)*cos(eff_rot);
    end
    
    %% RESULTS FILE
    
    %% INITIALIZE VARIABLES
    
    % CONTROL POINTS
    CPX = zeros([N,1]);
    CPY = zeros([N,1]);
    % INFLUENCE MATRICES

    % TANGENTIAL VELOCITIES
    V_t = zeros([N,1]);

    %% CALCULATIONS
        
    %% GEOMETRIC PARAMETERS OF THE AIRFOIL
    
    % UPDATE COORDINATES OF CONTROL POINTS
    for i = 1:N
        CPX(i) = 0.5*(x(i)+x(i+1));
        CPY(i) = 0.5*(y(i)+y(i+1));
    end
    
    % AIRFOIL PERIMETER
    L = 0;
    for i = 1:N
        L = L + sqrt((x(i+1)-x(i))^2+(y(i+1)-y(i))^2);
    end
    
    % PANEL ANGLE
    theta_t = zeros([N,1]);
    for i = 1:N
        theta_t(i) = atan2(y(i+1)-y(i),x(i+1)-x(i));
    end
    % TANGENTIAL VECTOR
    t = zeros([N,2]);
    for i = 1:N
        t(i,1) = cos(theta_t(i));
        t(i,2) = sin(theta_t(i));
    end
    % NORMAL PANEL ANGLE AND VECTOR
    theta_n = theta_t + pi/2;
    n = zeros([N,2]);
    for i = 1:N
        n(i,1) = cos(theta_n(i));
        n(i,2) = sin(theta_n(i));
    end
    
    %% INFLUENCE COEFFICIENTS
    
    % A PANELS INFLUENCE
    [AN,AT] = calc_AN_AT(x,y,N);
    % B PANELS INFLUENCE
    [BN,BT] = calc_BN_BT(x,y,N);
    
    %% SYSTEM OF ALGEBRAIC EQUATIONS
    
    A1 = AN;
    A2 = sum(BN,2); % be careful
    A3 = AT(1,:)+AT(N,:);
    A4 = sum(BT(1,:)+BT(N,:));
    A = [A1 A2;A3 A4];
    
    b1 = -V_mag*sin(AoA-theta_t(:));
    %% KUTTA CONDITION
    b2 = -V_mag*cos(AoA-theta_t(1))-V_mag*cos(AoA-theta_t(N));
    b = [b1; b2];
    
    X = A\b;
    
    q = X(1:N);
    TAU = X(N+1);
    
    for i = 1:N
        term1 = sum(AT(i,:)*q);
        term2 = TAU*sum(BT(i,:));
        term3 = V_mag*cos(AoA-theta_t(i));
        V_t(i) = term1 + term2 + term3;
    end
    
    %% PRESSURE COEFFICIENT

    CP = 1-(V_t.^2)/V_mag^2;
    
    %% X FORCE COEFFICIENTS
    CX = 0;
    for i = 1:N
        CX = CX + CP(i)*(y(i+1)-y(i));
    end
    
    %% Y FORCE COEFFICIENTS
    CY = 0;
    for i = 1:N
        CY = CY - CP(i)*(x(i+1)-x(i));
    end
    
    %% MOMENT COEFFICIENTS
    C_M = 0;
    for i = 1:N
        C_M = C_M + CP(i)*((x(i+1)-x(i))*CPX(i)+(y(i+1)-y(i))*CPY(i));
    end
    
    %% DRAG, LIFT AND MOMENT COEFFICIENTS
    CD = CX*cos(AoA)+CY*sin(AoA);
    CL = CY*cos(AoA)-CX*sin(AoA);
    CM = C_M;

end

%% PLOTS FOR CALCULATIONS

vtx = V_t.*cos(theta_t);
vty = V_t.*sin(theta_t);
cprx = -CP.*cos(theta_n);
cpry = -CP.*sin(theta_n);
quiver(CPX,CPY,vtx,vty);
hold on;
quiver(CPX,CPY,cprx,cpry);
daspect([1 1 1]);

%% END
WarnWave = [sin(1:0.25:500), sin(1:0.5:500), sin(1:1.0:500)];
Audio = audioplayer(WarnWave, 11025);
play(Audio);
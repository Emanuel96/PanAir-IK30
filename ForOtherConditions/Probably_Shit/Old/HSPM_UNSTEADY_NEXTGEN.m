%  UNSTEADY HSPM by Komasho (Based on TENG)
%  Last Update: June 17th, 2022 by Emanuel Camacho

clc; close all; %clear all;

%% AIRFOIL COORDINATES
% (COORDINATES MUST ALWAYS BE CREATED/READ ON THE COUNTERCLOCKWISE DIRECTION)

% AIRFOIL CHORD
chord = 0.2;
% NUMBER OF BOUNDARY POINTS (ODD NUMBER)
N_P = 31;
if rem(N_P,2) == 0 % FORCE ODD NUMBER
    N_P = N_P + 1;
end
% NUMBER OF PANELS
N = N_P-1;

%% NACA0012 COORDINATE GENERATOR

% COSSINE SPACING
beta = linspace(0,2*pi,N_P)';
x1 = 0.5*(1-cos(beta+pi));
y1 = 0.594689181*(0.298222773*sqrt(x1)-0.127125232.*x1-0.357907906*x1.^2+0.291984971*x1.^3-0.105174606*x1.^4);
pivot = 0.3;
x = (x1-pivot)*chord;
y=([y1(1:(N_P-1)/2);-y1((N_P-1)/2+1:N_P)])*chord;
y=flip(y);
x=round(x,7);
y=round(y,7);

% AIFOIL COORDINATES WITHOUT DEFORMATION
x_wd = x;
y_wd = y;

% CHECK LEADING EDGE COORDINATES
for i = 1:N
    if x(i) > 0 && x(i+1) < 0
        iSLE = i;
    elseif x(i) < 0 && x(i+1) > 0
        iELE = i+1;
    end
end

% UNCOMMENT IF THE WHOLE AIRFOIL IS PITCHING
%iSLE = -1000;
%iELE = 1000;

%% SIMULATION PARAMETERS

% GOVERNING PARAMETERS
AoA = 0*pi/180;
Re = 1e4;
NonDimAmpl = 0.25;
RedFreq = 4;
%AMP_AL = -atan(RedFreq*NonDimAmpl);
AMP_AL = 0*pi/180;
% NUMBER OF SIMULATED PERIODS
NPs = 5;
% NUMBER OF TIME STEPS PER PERIOD
NTSPP = 100;

zeta = 0.5;

%% SIMULATION PARAMETERS (CONTINUATION)

% VELOCITY CALCULATED BASED ON THE REYNOLDS NUMBER
V_INF = [Re*1.7894e-5/(chord*1.225)*cos(AoA) Re*1.7894e-5/(chord*1.225)*sin(AoA)];
V_mag = (V_INF(1)^2+V_INF(2)^2)^0.5;

% AMPLITUDE AND FREQUENCY
AMP = NonDimAmpl*chord;
freq = RedFreq*V_mag/(2*pi*chord);
T = 1/freq;

% TIME STEP CALCULATION (IN SECONDS)
dt = (1/freq)/NTSPP;
% NUMBER OF TIMESTEPS
NTS = NPs*NTSPP;
% MAX NUMBER OF ITERATIONS PER TIMESTEP
NITER = 50;
% MAX NUMBER OF VORTICES
NVs = 2*NTSPP;

%% RESULTS FILE

fid = fopen('Results.dat','w');
fprintf(fid,'VARIABLES=X,Y,VX,VY,CV\n');

fvec = fopen('PressAndVel.dat','w');
fprintf(fvec,'VARIABLES=X,Y,VX,VY,CPX,CPY,CP\n');

f_res = fopen('residuals.dat','w');
cont = 1;

%% OPTIMIZATION FILES

fop = fopen('Results_Opt.dat','w');
fprintf(fop,'VARIABLES=i,AMP,GAMMA\n');

AMP_AL_old = AMP_AL/2;
AMP_AL_old_1 = AMP_AL/4;
CP_old = 0;
CP_old_1 = 0;
gamma = 0.01;

%% INITIALIZE VARIABLES

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
AN = zeros([N,N]);
AT = zeros([N,N]);
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

%% VARIABLES FOR PRESSURE DISTRIBUTION

%EXPANSION RATIO (GEOMETRIC SERIES)
RATIO = 1.1;
%NUMBER OF PANELS IN INTEGRATION PATH
N_PHI = 20;

XZ = zeros([N_PHI+1,1]);
YZ = zeros([N_PHI+1,1]);
XMIDZ = zeros([N_PHI,1]);
YMIDZ = zeros([N_PHI,1]);
V_T_PHI_F = zeros([N_PHI,1]);
V_T_PHI = zeros([N,1]);
PHI_NODE = zeros([N+1,1]);
PHI = zeros([N,1]);
CP = zeros([N,1]);

%% FORCES

CD = zeros([NTS,1]);
CL = zeros([NTS,1]);
CM = zeros([NTS,1]);

V_t = zeros([N,1]);

V_XY_AIRFOIL = zeros([N,2]);
V_AL_AIRFOIL = zeros([N,2]);
V_AL_NODE_AIRFOIL = zeros([N+1,2]);

V_STREAM = zeros([N,2]);

b = zeros([N,1]);
c = zeros([N,1]);

q = zeros([N,1]);

%% INITIALIZE TAU, THETA AND DELTA

UW_old = V_INF(1);
VW_old = V_INF(2);
UW = V_INF(1);
VW = V_INF(2);
TAUM1 = 0;
TAU = 0;
TAU_old = 0;

HXO = 0;
HYO = 0;

PHIM1 = zeros([N,1]);

Theta = atan2(VW,UW);
Delta = dt*sqrt(UW^2+VW^2);

%WAKE ELEMENT
x(N+2) = x(N+1)+Delta*cos(Theta);
y(N+2) = y(N+1)+Delta*sin(Theta);

NM = 0; %%NUMBER OF CORE VORTICES

%% CALCULATIONS

iter = 0;

for k = 1:NTS
    
    %% COMPUTE RIGID BODY MOTION
       
    for i = 1:N
        
        % PLUNGING
%         dt = (1/freq)/(NTSPP);
%         HY = AMP*cos(2*pi*freq*dt*k)-AMP;
%         DHY = HY-HYO;
%         V_XY_AIRFOIL(i,:) = [0 -2*pi*freq*AMP*sin(2*pi*freq*dt*k)];
%         PLUNG_VEL(k) = -V_XY_AIRFOIL(i,2);
        
        if k >= floor(k/NTSPP)*NTSPP && k < floor(k/NTSPP)*NTSPP+NTSPP/2
            dt = zeta*(1/freq)/(NTSPP/2);
            HY = AMP*cos(pi*freq/zeta*(dt*k))-AMP;
            V_XY_AIRFOIL(i,:) = [0 -pi*freq*AMP/zeta*sin(pi*freq/zeta*(dt*k))];
        elseif k >= floor(k/NTSPP)*NTSPP+NTSPP/2 && k < floor(k/NTSPP)*NTSPP+NTSPP
            dt = (1-zeta)*(1/freq)/(NTSPP/2);
            HY = AMP*cos(pi*freq/(1-zeta)*(T-dt*k))-AMP;
            V_XY_AIRFOIL(i,:) = [0 pi*freq*AMP/(1-zeta)*sin(pi*freq/(1-zeta)*(T-dt*k))];
        end
        DHY = HY-HYO;
        PLUNG_VEL(k) = -V_XY_AIRFOIL(i,2);
                
        % PITCHING
        V_AL_NODE_AIRFOIL(i,:) = [0 0];
        V_AL_AIRFOIL(i,:) = [0 0];
        ROT(k) = 0;
        
        if i > iSLE && i < iELE 
            OMEGA = -[0 0 -2*pi*freq*AMP_AL*sin(2*pi*freq*dt*k+pi/2)];
            %OMEGA = -[0 0 -2*pi*freq*AMP_AL*sin(2*pi*freq*dt*k+pi/2)];
            ROT(k) = OMEGA(3);

            %AT AIRFOIL NODES
            vecpos_node = [x(i) y(i) 0];
            vel_node = cross(OMEGA,vecpos_node);
            V_AL_NODE_AIRFOIL(i,:) = [vel_node(1) vel_node(2)];

            %AT THE CONTROL POINTS
            vecpos = [0.5*(x(i)+x(i+1)) 0.5*(y(i)+y(i+1)) 0];
            vel = cross(OMEGA,vecpos);
            V_AL_AIRFOIL(i,:) = [vel(1) vel(2)]; 
        end
    end
    
    % DEFORMATION IN BODY AXIS
    for i = 1:N+1
        if i > iSLE && i < iELE 
            x(i) = x_wd(i)*cos(AMP_AL*sin(2*pi*freq*dt*k))-y_wd(i)*sin(AMP_AL*sin(2*pi*freq*dt*k)); 
            y(i) = x_wd(i)*sin(AMP_AL*sin(2*pi*freq*dt*k))+y_wd(i)*cos(AMP_AL*sin(2*pi*freq*dt*k));
        end  
    end

    % COORDINATES OF CONTROL POINTS
    for i = 1:N
        CPX(i) = 0.5*(x(i)+x(i+1));
        CPY(i) = 0.5*(y(i)+y(i+1));
    end
    
    %% PANEL TANGENTIAL AND NORMAL ANGLES
    
    theta_t = zeros([N,1]);
    for i = 1:N
         theta_t(i) = atan2(y(i+1)-y(i),x(i+1)-x(i));
    end
    
    t = zeros([N,2]);
    for i = 1:N
        t(i,1) = cos(theta_t(i)); 
        t(i,2) = sin(theta_t(i));
    end
    
    theta_n = theta_t + pi/2;
    n = zeros([N,2]);
    for i = 1:N
        n(i,1) = cos(theta_n(i)); 
        n(i,2) = sin(theta_n(i));
    end
        
    %% VELOCITY STREAM

    for i = 1:N
        V_STREAM(i,:) = V_INF + V_XY_AIRFOIL(i,:);
    end
       
    %% AIRFOIL PERIMETER
    
    L = 0;
    for i = 1:N
        L = L + sqrt((x(i+1)-x(i))^2+(y(i+1)-y(i))^2);
    end

    %% CONVECT CORE VORTICES DOWNSTREAM

    if NM >= 1
        for m = 1:NM
            XM(m) = XM(m) + dt*UM(m);
            YM(m) = YM(m) + dt*VM(m) + DHY;
        end
    end
    
    %% ITERATIVE CALCULATIONS TO OBTAIN DELTA, THETA AND TAU_W
    
    for iter = 1:NITER

        %% DISTANCE BETWEEN VORTICES
        % Find all neighbors within specified distance using input data
        %[Idx,D] = rangesearch([x y],[x y],r);
        
        %% A INFLUENCE COEFFICIENTS
      
        %A PANELS INFLUENCE
        [AN,AT] = calc_AN_AT(x,y,N);
        %A SHED VORTICITY INFLUENCE
        [AYW_J,AXW_J] = calc_AXW_J_AYW_J(x,y,N);
        %A CORE VORTEX INFLUENCE
        if NM >= 1
            [AYH_J,AXH_J] = calc_AXH_J_AYH_J(x,y,XM,YM,N,NM);
        end
     
        %% B INFLUENCE COEFFICIENTS

        %B PANELS INFLUENCE
        [BN,BT] = calc_BN_BT(x,y,N);
        %B SHED VORTICITY INFLUENCE
        [BNI_W,BTI_W] = calc_BNI_W_BTI_W(x,y,N);
        [BYW_J,BXW_J] = calc_BXW_J_BYW_J(x,y,N);
        if NM >= 1
            [BYH_W,BXH_W] = calc_BXH_W_BYH_W(x,y,XM,YM,N,NM);
        end
        %B CORE VORTEX INFLUENCE
        if NM >= 1
            [BYH_J,BXH_J] = calc_BXH_J_BYH_J(x,y,XM,YM,N,NM);
        end

        %% C INFLUENCE COEFFICIENTS

        if NM >= 1
            %CORE VORTEX INFLUENCE ON PANELS
            [CNI_M,CTI_M] = calc_CNI_M_CTI_M(x,y,XM,YM,N,NM);
            %CORE VORTEX INFLUENCE ON WAKE
            [CYW_M,CXW_M] = calc_CXW_M_CYW_M(x,y,XM,YM,N,NM);
            %CORE CORTEX INFLUENCE ON OTHER VORTEX
            [CYH_M,CXH_M] = calc_CXH_M_CYH_M(XM,YM,NM);
        end
        
        %% SYSTEM OF ALGEBRAIC EQUATIONS
        
        %A MATRIX
        A = AN;

        for i = 1:N
            %B ARRAY
            term1 = L/Delta*BNI_W(i);
            term2 = -sum(BN(i,:));
            b(i) = term1+term2;
            % C ARRAY
            term1 = -(dot(V_STREAM(i,:),n(i,:)));
            term2 = -L/Delta*TAUM1*BNI_W(i);
            term3 = 0; 
            if NM >= 1
                for m = 1:NM
                    term3 = term3 - CNI_M(i,m)*CV(m);
                end
            end
            c(i) = term1 + term2 + term3 + dot(V_AL_AIRFOIL(i,:),n(i,:)); % NEW FORMULATION
        end

        %% SOLVE LINEAR SYSTEM AS A FUNCTION OF VORTICITY DISTRIBUTION
        
        % A*q = TAU*B + C;
        % q = b1*TAU + b2;
        b1 = A\b;
        b2 = A\c;
        
        %% CALCULATE TANGENTIAL VELOCITIES AS A FUNCTION OF VORTICITY DISTRIBUTION
        
        % D11 TERM
        term1 = sum(AT(1,:).*b1')+sum(BT(1,:));
        term2 = -L/Delta*BTI_W(1);
        D11 = term1 + term2;
        % D21 TERM
        term1 = sum(AT(1,:).*b2');
        term2 = L*TAUM1/Delta*BTI_W(1);
        term3 = dot(V_STREAM(1,:),t(1,:)); %CHECK THIS THING
        term4 = 0;
        if NM >= 1
            for m = 1:NM
                term4 = term4 + CTI_M(1,m)*CV(m);
            end
        end
        D21 = term1 + term2 + term3 + term4;
        
        % D1N TERM
        term1 = sum(AT(N,:).*b1')+sum(BT(N,:));
        term2 = -L/Delta*BTI_W(N);
        D1N = term1 + term2;
        % D2N TERM
        term1 = sum(AT(N,:).*b2');
        term2 = L*TAUM1/Delta*BTI_W(N);
        term3 = dot(V_STREAM(N,:),t(N,:)); % CHECK THIS THING
        term4 = 0;
        if NM >= 1
            for m = 1:NM
                term4 = term4 + CTI_M(N,m)*CV(m);
            end
        end
        D2N = term1 + term2 + term3 + term4;
        
        %% KUTTA CONDITION 
        
        % VERIFY THESE TERMS
        C1 = D11^2-D1N^2;
        C2 = 2*D11*D21-2*D1N*D2N-2*L/dt;
        C3 = D21^2-D2N^2+2*L*TAUM1/dt;
               
        RADI = sqrt(C2^2-4*C1*C3);
        
        TAU = (-C2-RADI)/(2*C1);
                        
        %% UPDATED Q VECTOR WITH NEW TAU FROM KUTTA CONDITION
        q = TAU*(A\b) + A\c;
        
        %% TANGENTIAL VELOCITIES
        
        for i = 1:N
            term1 = sum(AT(i,:).*q');
            term2 = TAU*sum(BT(i,:));
            term3 = dot(V_STREAM(i,:),t(i,:));
            term4 = L*(TAUM1-TAU)/Delta*BTI_W(i);
            term5 = 0;
            if NM >= 1
                for m = 1:NM
                    term5 = term5 + CTI_M(i,m)*CV(m);
                end
            end
            V_t(i) = term1 + term2 + term3 + term4 + term5;
        end  
        
        %% UPDATE UW AND VW

        [UW,VW] = UW_VW(AXW_J,AYW_J,q,BXW_J,BYW_J,TAU,V_INF,NM,CXW_M,CYW_M,CV,k);
                
        %% NEW VALUES FOR NEXT ITERATION AND CHECK CONVERGENCE
        
        Theta = atan2(VW,UW);
        Delta = dt*sqrt(UW^2+VW^2);
                
        x(N+2) = x(N+1)+Delta*cos(Theta);
        y(N+2) = y(N+1)+Delta*sin(Theta);
        
        %fprintf('Timestep: %3d | Iteration: %2d (Delta: %2.3fm | Theta: %2.5fº)\n',k,iter,Delta,Theta*180/pi);
        
        if abs(VW - VW_old)/abs(VW_old) <= 0.01 && abs(UW - UW_old)/abs(UW_old) <= 0.01
            break;
        end
        
        cont = cont + 1;
        fprintf(f_res,'%2.5f %2.5f %2.5f\n',cont,abs(VW - VW_old),abs(VW - VW_old));
        
        UW_old = UW;
        VW_old = VW;
        
        if iter == NITER
            fprintf('Convergence not achieved. Change time step and/or number of panels.');
            return;
        end
        
    end
    
    %% COMPUTE AERO COEFFICIENTS
    
    %i_LE = N_P/2+1;
    i_LE = (N_P-1)/2+1;
    XZ(1) = x(i_LE);
    YZ(1) = y(i_LE);
    
    % TEN CHORDS AWAY FROM LEADING EDGE
    AX = (10*(1-RATIO))/(1-RATIO^(N_PHI+1));
    
    % ANGLE OF ATTACK
    alpha = atan(V_INF(2)/V_INF(1));
    
    for z = 2:N_PHI+1
        XZ(z) = XZ(z-1)-(AX*RATIO^(z-1))*cos(alpha);
        YZ(z) = YZ(z-1)-(AX*RATIO^(z-1))*sin(alpha);
    end
    
    for z = 1:N_PHI
        XMIDZ(z) = 0.5*(XZ(z)+XZ(z+1));
        YMIDZ(z) = 0.5*(YZ(z)+YZ(z+1));
    end
    
    %% TANGENTIAL VELOCITIES AT Z PANELS
    
    [AT_F,BT_F,BTF_W] = MATRICES_FOR_FORCES(N_PHI,N,XZ,YZ,x,y,alpha);
    
    if NM >= 1
        CTF_M = MATRICES_FOR_FORCES_C(N_PHI,XZ,YZ,alpha,XM,YM,NM);
    end
   
    for f = 1:N_PHI
        term1 = sum(AT_F(f,:).*q');
        term2 = TAU*sum(BT_F(f,:));
        term3 = L*(TAUM1-TAU)/Delta*BTF_W(f);
        term4 = 0;
        if NM >= 1
            for m = 1:NM
                term4 = term4 + CTF_M(f,m)*CV(m);
            end
        end
        V_T_PHI_F(f) = term1 + term2 + term3 + term4;
    end  
    
    %% PHI AT THE LEADING EDGE
    PHI_LE = 0;
    for f = 1:N_PHI
        PHI_LE = PHI_LE + V_T_PHI_F(f)*((XZ(f+1)-XZ(f))^2+(YZ(f+1)-YZ(f))^2)^0.5;
    end
    
    %% TANGENTIAL VELOCITIES FOR PHI CALCULATION
    for i = 1:N
        term1 = sum(AT(i,:).*q');
        term2 = TAU*sum(BT(i,:));
        term3 = L*(TAUM1-TAU)/Delta*BTI_W(i);
        term4 = 0;
        if NM >= 1
            for m = 1:NM
                term4 = term4 + CTI_M(i,m)*CV(m);
            end
        end
        V_T_PHI(i) = term1 + term2 + term3 + term4;
    end  
    
    %LOWER SURFACE
    for i = 1:i_LE-1
        term1 = 0;
        for j = i:i_LE-1
            term1 = term1 + V_T_PHI(j)*((x(j+1)-x(j))^2+(y(j+1)-y(j))^2)^0.5;
        end
        PHI_NODE(i) = PHI_LE - term1;
    end
    
    %UPPER SURFACE
    for i = i_LE:N+1
        term1 = 0;
        for j = i_LE:i-1
            term1 = term1 + V_T_PHI(j)*((x(j+1)-x(j))^2+(y(j+1)-y(j))^2)^0.5;
        end
        PHI_NODE(i) = PHI_LE + term1;
    end
    
    %% VELOCITY POTENTIAL
    for i = 1:N
        PHI(i) = 0.5*(PHI_NODE(i)+PHI_NODE(i+1));
    end
    
    %% PRESSURE COEFFICIENT
    for i = 1:N
        term1 = (V_STREAM(i,1)^2+V_STREAM(i,2)^2)/V_mag^2;
        term2 = (V_t(i)^2+dot(V_AL_AIRFOIL(i,:),n(i,:))^2)/V_mag^2; % NEW FORMULATION
        term3 = (2/V_mag^2)*(PHI(i)-PHIM1(i))/dt;
        CP(i) = term1 - term2 - term3;
    end
   
    %% X FORCE COEFFICIENTS
    CX = 0;
    for i = 1:N
       CX = CX + CP(i)*(y(i+1)-y(i));
       %panel_length = sqrt((x(i+1)-x(i))^2+(y(i+1)-y(i))^2);
       %CX = CX + CP(i)*(y(i+1)-y(i)) + sign(V_t(i)*cos(theta_t(i)))*1.328/(panel_length*1.225*abs(V_t(i))/1.7894e-5)^0.5*panel_length*cos(theta_t(i));
    end

    %% Y FORCE COEFFICIENTS
    CY = 0;
    for i = 1:N
       CY = CY - CP(i)*(x(i+1)-x(i));
    end

    %% MOMENT COEFFICIENTS
    C_M = 0;
    for i = 1:N
        if i > iSLE && i < iELE
            C_M = C_M + CP(i)*((x(i+1)-x(i))*CPX(i)+(y(i+1)-y(i))*CPY(i));
        end
    end

    %% DRAG, LIFT AND MOMENT COEFFICIENTS
    CD(k) = CX*cos(alpha)+CY*sin(alpha);
    CL(k) = CY*cos(alpha)-CX*sin(alpha);
    CM(k) = C_M;

    %% VORTEX CORE VELOCITIES
    
    if NM >= 1
        [UM,VM] = UM_VM(AXH_J,AYH_J,q,BXH_J,BYH_J,BXH_W,BYH_W,TAU,TAUM1,Delta,L,V_INF,NM,CXH_M,CYH_M,CV,k);
    end
    
    %% DEATTACH TRAILING-EDGE PANEL AND CREATE CORE VORTEX
    
    NM = NM + 1;
    
    if NM == NVs
        NM = NVs-1;
        UM = circshift(UM,-1);
        VM = circshift(VM,-1);
        XM = circshift(XM,-1);
        YM = circshift(YM,-1);
        CV = circshift(CV,-1);
    end
       
    UM(NM) = UW;
    VM(NM) = VW;
    XM(NM) = x(N+1) + 0.5*Delta*cos(Theta);
    YM(NM) = y(N+1) + 0.5*Delta*sin(Theta);
    % CORE VORTEX CIRCULATION
    CV(NM) = L*(TAUM1-TAU);

    %% PREPARATIONS FOR THE NEXT TIMESTEP
   
    HYO = HY;
    TAUM1 = TAU;
    PHIM1 = PHI;

    %% PLOTS FOR CALCULATIONS

    vtx = V_t.*cos(theta_t);
    vty = V_t.*sin(theta_t);
    cprx = -CP.*cos(theta_n);
    cpry = -CP.*sin(theta_n);

    %% DATA FOR TECPLOT
    
    % PRESURE AND VELOCITY DISTRIBUTION
    
    fprintf(fvec,"ZONE T='RESULTS' I=%d\n",N);

    for airfoil_c = 1:N
        fprintf(fvec,"%f %f %f %f %f %f %f\n",CPX(airfoil_c),CPY(airfoil_c),vtx(airfoil_c),vty(airfoil_c),cprx(airfoil_c),cpry(airfoil_c),CP(airfoil_c));
    end
    
    % VORTEX
    
    fprintf(fid,"ZONE T='RESULTS' I=%d\n",NM+N+1);
    
    for airfoil_c = 1:N+1
        fprintf(fid,"%f %f %f %f %f\n",x(airfoil_c),y(airfoil_c),V_AL_NODE_AIRFOIL(airfoil_c,1),V_AL_NODE_AIRFOIL(airfoil_c,2),0);
    end
    
    for vor = 1:NM
        fprintf(fid,"%f %f %f %f %f\n",XM(vor),YM(vor),UM(vor),VM(vor),CV(vor));
    end
    
%     pause(0.1);
%     subplot(2,2,1)
%     plot(atan(PLUNG_VEL(1:k)/V_mag),CT(1:k));
%     subplot(2,2,2)
%     plot(atan(PLUNG_VEL(1:k)/V_mag),CL(1:k));
%     subplot(2,2,3)
%     plot(atan(PLUNG_VEL(1:k)/V_mag),CM(1:k));
%     subplot(2,2,4)
%     plot(atan(PLUNG_VEL(1:k)/V_mag),CPR(1:k));
    
    %OPTIMIZATION - LEADING-EDGE PITCHING AMPLITUDE
    
%     if rem(k+1,NTSPP) == 0
%         %THRUST COEFFICIENT
%         CT = -CD;
%         %POWER COEFFICIENT
%         plungingP = -PLUNG_VEL.*CL/V_mag;
%         CM = -CM; %VERIFY THIS
%         pitchingP = -ROT.*CM*chord/V_mag;
%         CPR = plungingP+pitchingP;
%         %PROPULSIVE EFFICIENCY
%         %ETA = mean(CT(NTS-NTSPP:NTS))/mean(CPR(NTS-NTSPP:NTS));
%     
%         CP = mean(CPR(k-NTSPP+2:k));
%                
%         AMP_AL = AMP_AL_old - gamma*(CP_old-CP_old_1)/(AMP_AL_old-AMP_AL_old_1);
%         
%         fprintf(fop,'%f %f\n',CP,AMP_AL*180/pi);
%         
%         AMP_AL_old_1 = AMP_AL_old;
%         AMP_AL_old = AMP_AL;
%         CP_old_1 = CP_old;
%         CP_old = CP;
%     end
    
    if rem(k,NTSPP) == 0
        %% MEAN DRAG
        meanCT = mean(-CD(k-NTSPP+1:k));
        %% MEAN POWER
        plungingP = -PLUNG_VEL.*CL/V_mag;
        meanCP = mean(plungingP(k-NTSPP+1:k));
        %% PROPULSIVE EFFICIENCY
        ETA = meanCT/meanCP;
        
        fprintf('Timestep: %3d | meanCT: %2.3f | meanCP: %2.3f | Theta: %2.3f\n',k,meanCT,meanCP,ETA);
    end

end

fclose(fid);
fclose(fvec);
fclose(fop);

fclose(f_res);

%% PROPULSIVE CALCULATIONS

%THRUST COEFFICIENT (WITH VISCOUS CORRECTION)
%CT = -(CD+0.455/log10(Re)^2.58*L/chord);
%THRUST COEFFICIENT (WITH VISCOUS CORRECTION)
CT = -CD;
%POWER COEFFICIENT
plungingP = -PLUNG_VEL.*CL/V_mag;
CM = -CM; %VERIFY THIS
pitchingP = -ROT.*CM*chord/V_mag;
CPR = plungingP+pitchingP;
%PROPULSIVE EFFICIENCY
ETA = mean(CT(NTS-NTSPP:NTS))/mean(CPR(NTS-NTSPP:NTS));

%% FORCES DISPLAY

fprintf('\n======= RESULTS =======\n');
fprintf('Mean Thrust Coefficient (CT) : %2.3f\n',mean(CT(NTS-NTSPP:NTS)));      
fprintf('Mean Power Coefficient (CPR) : %2.3f\n',mean(CPR(NTS-NTSPP:NTS)));  
fprintf('Propulsive Efficiency (ETA) : %2.3f\n',ETA);

%% FORCES FILE

fid_force = fopen('forces.dat','w');
%fprintf(fid_force,'VARIABLES=t/T,VEL,AoA<sub>eff</sub>,CD,CL,CM,CT,CPR\n');
fprintf(fid_force,'t y_dot alpha_eff cd cl cm ct cpr\n');

for k = 2:NTS
    fprintf(fid_force,"%f %f %f %f %f %f %f %f\n",k*dt/(1/freq),PLUNG_VEL(k),atan(PLUNG_VEL(k)/V_mag)*180/pi,CD(k),CL(k),CM(k),CT(k),CPR(k));
end

% subplot(3,1,1)
% plot(CT);
% subplot(3,1,2)
% plot(CL);
% subplot(3,1,3)
% plot(CM);

figure
plot(x,y);
axis([-1 11 -3 3]);
pbaspect([2 1 1]);
hold on
scatter(XM,YM,5,CV,'filled');

figure
subplot(2,2,1)
plot(atan(PLUNG_VEL(:)/V_mag)*180/pi,CT);
subplot(2,2,2)
plot(atan(PLUNG_VEL(:)/V_mag)*180/pi,CL);
subplot(2,2,3)
plot(atan(PLUNG_VEL(:)/V_mag)*180/pi,CM);
subplot(2,2,4)
plot(atan(PLUNG_VEL(:)/V_mag)*180/pi,CPR);

fclose(fid_force);

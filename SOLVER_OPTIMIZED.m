%  UNSTEADY HSPM by Komasho (Based on TENG)
%  Last Update: MARCH 28th, 2022 by Emanuel Camacho

tic

clc; close all; %clear all;

%% SIMULATION PARAMETERS
% AERODYNAMIC CHORD
chord = 1;
% ANGLE OF ATTACK
AoA = 0*pi/180;
% REYNOLDS NUMBER
Re = 1e4;
%Re = [1e4]';
% NONDIMENSIONAL AMPLITUDE
NonDimAmpl = 0.5;
%NonDimAmpl = [0.1 0.2 0.4]';
% REDUCED FREQUENCY
RedFreq = 0.1;
%RedFreq = [0.25:0.25:0.75 1:0.5:4]';
%RedFreq = [2:2:8 10:2.5:15]';
% MEAN ANGLE OF ATTACK
alpha_0 = -0*pi/180;
alpha_1 = -0*pi/180;
% PITCHING AMPLITUDE
%PitchAmp = -atand(NonDimAmpl*RedFreq)*pi/180;
PitchAmp = -0*pi/180;
PitchAmp_2 = -0*pi/180;
%PitchAmp = -[0 10 20]'*pi/180;
%PitchAmp_2 = -[0]'*pi/180;
% PIVOT (CHORD PERCENTAGE)
pivot = 0.3;
% ONLY LEADING EDGE?
onlyLE = 0; % 0 - STANDARD FLAPPING | 1 - LEADING-EDHE | 2 - FLAP
% ASYMMETRY OF KINEMATICS
zeta = 0.5;
%zeta = [0.5]';

for iRe = 1:size(Re)
    for iNonDimAmpl = 1:size(NonDimAmpl)
        for iRedFreq = 1:size(RedFreq)
            for iA_alpha = 1:size(PitchAmp)
                for iA_alpha_2 = 1:size(PitchAmp_2)
                    for iZeta = 1:size(zeta)
                        
                        % AIRFOIL VELOCITY AND MAGNITUDE
                        V_INF = [Re(iRe)*1.7894e-5/(chord*1.225)*cos(AoA) Re(iRe)*1.7894e-5/(chord*1.225)*sin(AoA)];
                        V_mag = (V_INF(1)^2+V_INF(2)^2)^0.5;
                        % MOTION AMPLITUDE
                        AMP = NonDimAmpl(iNonDimAmpl)*chord;
                        % PITCHING AMPLITUDE
                        AMP_AL = PitchAmp(iA_alpha);
                        % PITCHING AMPLITUDE - MAIN PART
                        AMP_AL_2 = PitchAmp_2(iA_alpha_2);
                        % MOTION FREQUENCY
                        freq = RedFreq(iRedFreq)*V_mag/(2*pi*chord);
                        %freq = RedFreq(iRedFreq)*V_mag/(pi*chord);
                        % MOTION PERIOD
                        T = 1/freq;
                        
                        fprintf('======= CURRENT CONDITION =======\n');
                        fprintf('k: %2.3f | h: %2.3f | A: %2.3f | A_2: %2.3f\n',RedFreq(iRedFreq),NonDimAmpl(iNonDimAmpl),-PitchAmp(iA_alpha)*180/pi,-PitchAmp_2(iA_alpha_2)*180/pi);
                        
                        folder_name = sprintf('Results/Re_%6.0f_h%0.3f_k%0.2f_A%0.2f',Re(iRe),NonDimAmpl(iNonDimAmpl),RedFreq(iRedFreq),-PitchAmp(iA_alpha)*180/pi);
                        mkdir(folder_name);
                                                
                        %% NUMERICAL PARAMETERS
                        
                        % NUMBER OF SIMULATED PERIODS
                        NPs = 10;
                        % NUMBER OF TIME STEPS PER PERIOD
                        NTSPP = 100;
                        % TIME STEP CALCULATION (IN SECONDS)
                        dt = (1/freq)/NTSPP;
                        % NUMBER OF TIMESTEPS
                        NTS = NPs*NTSPP;
                        % MAX NUMBER OF ITERATIONS PER TIMESTEP
                        NITER = 100;
                        % MAX NUMBER OF VORTICES
                        NVs = 2*NTSPP;
                        % CONVERGENCE REPORT
                        reportConvergence = 5;
                        
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
                                                             
                        % AIFOIL COORDINATES WITHOUT DEFORMATION (+ROTATION)
                        x_wd = x;
                        y_wd = y;
                                                                        
                        % DETERMINE LEADING EDGE COORDINATES
                        if onlyLE == (1 || 2)
                            for i = 1:N
                                if x(i) > 0 && x(i+1) < 0
                                    iSLE = i;
                                elseif x(i) < 0 && x(i+1) > 0
                                    iELE = i+1;
                                end
                            end
                        else
                            iSLE = -1000;
                            iELE = 1000;
                        end
                                                
                        for i = 1:N+1
                            if i > iSLE && i < iELE
                                x_wd(i) = x(i)*cos(alpha_0)-y(i)*sin(alpha_0);
                                y_wd(i) = x(i)*sin(alpha_0)+y(i)*cos(alpha_0);
                            else
                                x_wd(i) = x(i)*cos(alpha_1)-y(i)*sin(alpha_1);
                                y_wd(i) = x(i)*sin(alpha_1)+y(i)*cos(alpha_1);
                            end
                        end
                                                
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
                        BN = zeros([N,N]);
                        BT = zeros([N,N]);
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
                        
                        %% FIND STEADY SOLUTION FIRST TO INITIALIZE SOURCE
                        % STRENGTH DISTRIBUTION AND VORTICITY
                        
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
                        
                        % A & B PANELS INFLUENCE
                        [AN,AT,BN,BT] = calc_AN_AT_BN_BT(x,y,N);
                        
                        %% SYSTEM OF ALGEBRAIC EQUATIONS
                        
                        A1 = AN;
                        A2 = sum(BN,2); % be careful
                        A3 = AT(1,:)+AT(N,:);
                        A4 = sum(BT(1,:)+BT(N,:));
                        SS1 = [A1 A2;A3 A4];
                        
                        b1 = -V_mag*sin(AoA-theta_t(:));
                        %% KUTTA CONDITION
                        b2 = -V_mag*cos(AoA-theta_t(1))-V_mag*cos(AoA-theta_t(N));
                        SS2 = [b1; b2];
                        
                        X = SS1\SS2;
                        
                        q = X(1:N);
                        TAU = X(N+1);
                        TAUM1 = TAU;
                        
                        for i = 1:N
                            term1 = sum(AT(i,:)*q);
                            term2 = TAU*sum(BT(i,:));
                            term3 = V_mag*cos(AoA-theta_t(i));
                            V_t(i) = term1 + term2 + term3;
                        end
                        
                        %% PRESSURE COEFFICIENT
                        CP = 1-(V_t.^2)/V_mag^2;
                        
                        for i = 1:N
                            CPX(i) = 0.5*(x(i)+x(i+1));
                            CPY(i) = 0.5*(y(i)+y(i+1));
                        end
                        vtx = V_t.*cos(theta_t);
                        vty = V_t.*sin(theta_t);
                        cprx = -CP.*cos(theta_n);
                        cpry = -CP.*sin(theta_n);
                                                
                        i_LE = (N_P-1)/2+1;
                        XZ(1) = x(i_LE);
                        YZ(1) = y(i_LE);
                        
                        % TEN CHORDS AWAY FROM LEADING EDGE
                        % EXPANSION RATIO (GEOMETRIC SERIES)
                        RATIO = 1.1;
                        AX = (10*(1-RATIO))/(1-RATIO^(N_PHI+1));
                            
                        for z = 2:N_PHI+1
                            XZ(z) = XZ(z-1)-(AX*RATIO^(z-1))*cos(AoA);
                            YZ(z) = YZ(z-1)-(AX*RATIO^(z-1))*sin(AoA);
                        end
                            
                        for z = 1:N_PHI
                            XMIDZ(z) = 0.5*(XZ(z)+XZ(z+1));
                            YMIDZ(z) = 0.5*(YZ(z)+YZ(z+1));
                        end
                        
                        %% TANGENTIAL VELOCITIES AT Z PANELS
                        
                        [AT_F,BT_F,BTF_W] = MATRICES_FOR_FORCES(N_PHI,N,XZ,YZ,x,y,AoA);
                        
                        for f = 1:N_PHI
                            term1 = sum(AT_F(f,:)*q);
                            term2 = TAU*sum(BT_F(f,:));
                            V_T_PHI_F(f) = term1 + term2;
                        end
                        
                        %% PHI AT THE LEADING EDGE
                        PHI_LE = 0;
                        for f = 1:N_PHI
                            PHI_LE = PHI_LE + V_T_PHI_F(f)*((XZ(f+1)-XZ(f))^2+(YZ(f+1)-YZ(f))^2)^0.5;
                        end
                        
                        %% TANGENTIAL VELOCITIES FOR PHI CALCULATION
                        for i = 1:N
                            term1 = sum(AT(i,:)*q);
                            term2 = TAU*sum(BT(i,:));
                            V_T_PHI(i) = term1 + term2;
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
                        for i = i_LE:N
                            term1 = 0;
                            for j = i_LE:i-1
                                term1 = term1 + V_T_PHI(j)*((x(j+1)-x(j))^2+(y(j+1)-y(j))^2)^0.5;
                            end
                            PHI_NODE(i) = PHI_LE + term1;
                        end
                        
                        %% VELOCITY POTENTIAL
                        for i = 1:N
                            PHIM1(i) = 0.5*(PHI_NODE(i)+PHI_NODE(i+1));
                        end
                        
                        %% TRANSIENT CALCULATION
                        
                        for k = 1:NTS
                            
                            %% GEOMETRIC PARAMETERS OF THE AIRFOIL
                            
                            for i = 1:N
                                % ASYMMETRIC PLUNGING
                                
                                if k < 0
                                    HY = 0;
                                    V_XY_AIRFOIL(i,:) = [0 0];
                                else
                                    if k >= floor(k/NTSPP)*NTSPP && k < floor(k/NTSPP)*NTSPP+NTSPP/2
                                        dt = zeta*(1/freq)/(NTSPP/2);
                                        HY = AMP*cos(pi*freq/zeta*(dt*k))-AMP;
                                        V_XY_AIRFOIL(i,:) = [0 -pi*freq*AMP/zeta*sin(pi*freq/zeta*(dt*k))];
                                    elseif k >= floor(k/NTSPP)*NTSPP+NTSPP/2 && k < floor(k/NTSPP)*NTSPP+NTSPP
                                        dt = (1-zeta)*(1/freq)/(NTSPP/2);
                                        HY = AMP*cos(pi*freq/(1-zeta)*(dt*k))-AMP;
                                        V_XY_AIRFOIL(i,:) = [0 -pi*freq*AMP/(1-zeta)*sin(pi*freq/(1-zeta)*(dt*k))];
                                    end
                                end
                                
                                % Y-COORDINATE DEVIATION
                                DHY = HY-HYO;
                                % PLUNGING VELOCITY
                                PLUNG_VEL(k) = -V_XY_AIRFOIL(i,2);
                                
                                % MOVE LEADING EDGE
                                if i > iSLE && i < iELE
                                    OMEGA = [0 0 2*pi*freq*AMP_AL*cos(2*pi*freq*dt*k)];
                                    ROT(k) = -OMEGA(3); %IS THE SIGN RIGHT? DANGEROUS
                                    % ROTATIONAL VELOCITY AT AIRFOIL NODES
                                    vecpos_node = [x(i) y(i) 0];
                                    vel_node = cross(OMEGA,vecpos_node);
                                    V_AL_NODE_AIRFOIL(i,:) = [vel_node(1) vel_node(2)];
                                    % ROTATIONAL VELOCITY AT THE CONTROL POINTS
                                    vecpos = [0.5*(x(i)+x(i+1)) 0.5*(y(i)+y(i+1)) 0];
                                    vel = cross(OMEGA,vecpos);
                                    V_AL_AIRFOIL(i,:) = [vel(1) vel(2)];
                                else
                                    OMEGA = [0 0 2*pi*freq*AMP_AL_2*cos(2*pi*freq*dt*k)];
                                    ROT(k) = -OMEGA(3);
                                    % ROTATIONAL VELOCITY AT AIRFOIL NODES
                                    vecpos_node = [x(i) y(i) 0];
                                    vel_node = cross(OMEGA,vecpos_node);
                                    V_AL_NODE_AIRFOIL(i,:) = [vel_node(1) vel_node(2)];
                                    % ROTATIONAL VELOCITY AT THE CONTROL POINTS
                                    vecpos = [0.5*(x(i)+x(i+1)) 0.5*(y(i)+y(i+1)) 0];
                                    vel = cross(OMEGA,vecpos);
                                    V_AL_AIRFOIL(i,:) = [vel(1) vel(2)];
                                end
                            end
                            
                            % Effective Angle of Attack               
                            ALPHA_1 = -(AMP_AL*sin(2*pi*freq*dt*k)+alpha_0); %POSITIVE UP
                            ALPHA_2 = -(AMP_AL_2*sin(2*pi*freq*dt*k)+alpha_1); %POSITIVE UP
                            ALPHA_3 = atan((pivot*sin(ALPHA_1)+(1-pivot)*sin(ALPHA_2))/(pivot*cos(ALPHA_1)+(1-pivot)*cos(ALPHA_2))); %ERROR
                            ALPHA_EFF(k) = atan(-PLUNG_VEL(k)/V_mag) + ALPHA_3;
                            ALPHA_EFF(k) = ALPHA_EFF(k)*180/pi;
                            
                            % DEFLECTION IN BODY AXIS (ONLY PITCHING)
                            for i = 1:N+1
                                if onlyLE == 0
                                    x(i) = x_wd(i)*cos(AMP_AL*sin(2*pi*freq*dt*k))-y_wd(i)*sin(AMP_AL*sin(2*pi*freq*dt*k));
                                    y(i) = x_wd(i)*sin(AMP_AL*sin(2*pi*freq*dt*k))+y_wd(i)*cos(AMP_AL*sin(2*pi*freq*dt*k));
                                elseif onlyLE == 1
                                    AMP_AL_1 = 1/(1+exp(50*(x_wd(i))))*AMP_AL;
                                    x(i) = x_wd(i)*cos(AMP_AL_1*sin(2*pi*freq*dt*k))-y_wd(i)*sin(AMP_AL_1*sin(2*pi*freq*dt*k));
                                    y(i) = x_wd(i)*sin(AMP_AL_1*sin(2*pi*freq*dt*k))+y_wd(i)*cos(AMP_AL_1*sin(2*pi*freq*dt*k));
                                elseif onlyLE == 2
                                    AMP_AL_11 = 1/(1+exp(-50*(x_wd(i))))*AMP_AL_2;
                                    x(i) = x_wd(i)*cos(AMP_AL_11*sin(2*pi*freq*dt*k))-y_wd(i)*sin(AMP_AL_11*sin(2*pi*freq*dt*k));
                                    y(i) = x_wd(i)*sin(AMP_AL_11*sin(2*pi*freq*dt*k))+y_wd(i)*cos(AMP_AL_11*sin(2*pi*freq*dt*k));
                                end
                            end
                            
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
                            
                            %% VELOCITY STREAM
                            
                            for i = 1:N
                                V_STREAM(i,:) = V_INF + V_XY_AIRFOIL(i,:);
                            end
                            
                            %% CONVECT CORE VORTICES DOWNSTREAM
                            
                            if NM >= 1
                                for m = 1:NM
                                    XM(m) = XM(m) + dt*UM(m);
                                    YM(m) = YM(m) + dt*VM(m) + DHY;
                                end
                            end
                                                       
                            % A AND B PANELS INFLUENCE
                            [AN,AT,BN,BT] = calc_AN_AT_BN_BT(x,y,N);
                            % A AND B CORE VORTEX INFLUENCE ON PANEL
                            if NM >= 1
                                [AYH_J,AXH_J,BYH_J,BXH_J] = calc_AXH_J_AYH_J_BXH_J_BYH_J(x,y,XM,YM,N,NM,NTS);
                            end
                            % CALCULATE C COEFFICIENTS IF THERE ARE VORTICES
                            if NM >= 1
                                %CORE VORTEX INFLUENCE ON PANELS
                                [CNI_M,CTI_M] = calc_CNI_M_CTI_M(x,y,XM,YM,N,NM);
                                %CORE CORTEX INFLUENCE ON OTHER VORTEX
                                [CYH_M,CXH_M] = calc_CXH_M_CYH_M(XM,YM,NM);
                            end
                            
                            %% ITERATIVE PROCEDURE TO CALCULATE WAKE ELEMENT
                            
                            for iter = 1:NITER
                                
                                %% INFLUENCE COEFFICIENTS
                                
                                % A AND B SHED VORTICITY INFLUENCE
                                [AYW_J,AXW_J,BYW_J,BXW_J] = calc_AXW_J_AYW_J_BXW_J_BYW_J(x,y,N);
                                % B CORE VORTEX INFLUENCE ON WAKE
                                if NM >= 1
                                    [BYH_W,BXH_W] = calc_BXH_W_BYH_W(x,y,XM,YM,N,NM,NTS);
                                end
                                % B SHED VORTICITY INFLUENCE
                                [BNI_W,BTI_W] = calc_BNI_W_BTI_W(x,y,N);
                                if NM >= 1
                                    %CORE VORTEX INFLUENCE ON WAKE
                                    [CYW_M,CXW_M] = calc_CXW_M_CYW_M(x,y,XM,YM,N,NM);
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
                                    c(i) = term1 + term2 + term3 + dot(V_AL_AIRFOIL(i,:),n(i,:));
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
                                term3 = dot(V_STREAM(1,:),t(1,:));
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
                                term3 = dot(V_STREAM(N,:),t(N,:));
                                term4 = 0;
                                if NM >= 1
                                    for m = 1:NM
                                        term4 = term4 + CTI_M(N,m)*CV(m);
                                    end
                                end
                                D2N = term1 + term2 + term3 + term4;
                                
                                %% KUTTA CONDITION
                                
                                C1 = D11^2-D1N^2;
                                C2 = 2*D11*D21-2*D1N*D2N-2*L/dt;
                                C3 = D21^2-D2N^2+2*L*TAUM1/dt;
                                RADI = sqrt(C2^2-4*C1*C3);
                                
                                % AIRFOIL VORTICITY
                                
                                if C1 == 0
                                    fprintf('Kutta condition may be compromised. Check C1 parameter.\n');
                                    C1 = 1;
                                end
                                
                                TAU = (-C2-RADI)/(2*C1);
                                                             
                                %% UPDATED Q VECTOR WITH NEW TAU FROM KUTTA CONDITION
                                
                                q = TAU*(A\b) + A\c;
                                
                                %% TANGENTIAL VELOCITIES
                                
                                for i = 1:N
                                    term1 = sum(AT(i,:)*q);
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
                                                               
                                % CHECK CONVERGENCE (ABSOLUTE)
                                if abs(VW - VW_old)/abs(VW_old) <= 0.01 && abs(UW - UW_old)/abs(UW_old) <= 0.01
                                    break;
                                end
                                                                
                                Theta_old = Theta;
                                Delta_old = Delta;
                                UW_old = UW;
                                VW_old = VW;
                                
                                if iter == NITER
                                    fprintf('Convergence not achieved. Change time step and/or number of panels.');
                                    return;
                                end
                                
                            end
                            
                            %% COMPUTE AERO COEFFICIENTS
                            
                            i_LE = (N_P-1)/2+1;
                            XZ(1) = x(i_LE);
                            YZ(1) = y(i_LE);
                            
                            % TEN CHORDS AWAY FROM LEADING EDGE
                            % EXPANSION RATIO (GEOMETRIC SERIES)
                            RATIO = 1.1;
                            AX = (10*(1-RATIO))/(1-RATIO^(N_PHI+1));
                            
                            for z = 2:N_PHI+1
                                XZ(z) = XZ(z-1)-(AX*RATIO^(z-1))*cos(AoA);
                                YZ(z) = YZ(z-1)-(AX*RATIO^(z-1))*sin(AoA);
                            end
                            
                            for z = 1:N_PHI
                                XMIDZ(z) = 0.5*(XZ(z)+XZ(z+1));
                                YMIDZ(z) = 0.5*(YZ(z)+YZ(z+1));
                            end
                            
                            %% TANGENTIAL VELOCITIES AT Z PANELS
                            
                            [AT_F,BT_F,BTF_W] = MATRICES_FOR_FORCES(N_PHI,N,XZ,YZ,x,y,AoA);
                            
                            if NM >= 1
                                CTF_M = MATRICES_FOR_FORCES_C(N_PHI,XZ,YZ,AoA,XM,YM,NM);
                            end
                            
                            for f = 1:N_PHI
                                term1 = sum(AT_F(f,:)*q);
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
                                term1 = sum(AT(i,:)*q);
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
                            for i = i_LE:N
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
                                if onlyLE == (1 || 2)
                                    term2 = (V_t(i)^2+dot(V_AL_AIRFOIL(i,:),n(i,:))^2)/V_mag^2;
                                else
                                    term2 = (V_t(i)^2)/V_mag^2;
                                end
                                term3 = (2/V_mag^2)*(PHI(i)-PHIM1(i))/dt;
                                CP(i) = term1 - term2 - term3;
                            end
                            
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
                                if i > iSLE && i < iELE % CALCULATE THE PITCHING COEFFCIIENT OF MOVABLE PANELS
                                    C_M = C_M + CP(i)*((x(i+1)-x(i))*CPX(i)+(y(i+1)-y(i))*CPY(i));
                                end
                            end
                            
                            %% DRAG, LIFT AND MOMENT COEFFICIENTS
                            CD(k) = CX*cos(AoA)+CY*sin(AoA);
                            CL(k) = CY*cos(AoA)-CX*sin(AoA);
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
                                                                                    
                            if rem(k,NTSPP) == 0 && floor(k/NTSPP) > 1
                                CT = -CD;
                                %POWER COEFFICIENT
                                plungingP = -PLUNG_VEL.*CL/V_mag;
                                %CM = -CM; %VERIFY THIS
                                pitchingP = -ROT.*CM*chord/V_mag;
                                CPR = plungingP+pitchingP;
                                
                                %% MEAN DRAG
                                meanCPP = -mean(CD(k-NTSPP:k));
                                %% MEAN LIFT
                                meanCPR = mean(CPR(k-NTSPP:k));
                                %% EFFICIENCY
                                ETA = meanCPP/meanCPR;
                                
                                fprintf('Timestep: %3d | meanCPP: %2.5f | meanCPR: %2.5f | ETA: %2.5f\n',k,meanCPP,meanCPR,ETA);
                            end
                            
                        end
                        
                        %% PROPULSIVE CALCULATIONS

                        %THRUST COEFFICIENT
                        CT = -CD;
                        %POWER COEFFICIENT
                        plungingP = -PLUNG_VEL.*CL/V_mag;
                        %CM = -CM; %VERIFY THIS
                        pitchingP = -ROT.*CM*chord/V_mag;
                        CPR = plungingP+pitchingP;
                        %PROPULSIVE EFFICIENCY
                        ETA = mean(CT(NTS-NTSPP:NTS))/mean(CPR(NTS-NTSPP:NTS));
                                                
                        %% FORCES DISPLAY
                        
                        fprintf('======= RESULTS =======\n');
                        fprintf('Mean Thrust Coefficient (CT) : %2.5f\n',mean(CT(NTS-NTSPP:NTS)));
                        fprintf('Mean Power Coefficient (CPR) : %2.5f\n',mean(CPR(NTS-NTSPP:NTS)));
                        fprintf('Propulsive Efficiency (ETA) : %2.5f\n\n',ETA);
                                                
                        %fprintf(fparstu,'%2.5f %2.5f %2.5f %2.5f %2.5f %2.5f %2.5f %2.5f\n',Re(iRe),NonDimAmpl(iNonDimAmpl),RedFreq(iRedFreq),PitchAmp(iA_alpha)*180/pi,PitchAmp_2(iA_alpha_2)*180/pi,mean(CT(NTS-NTSPP:NTS)),mean(CPR(NTS-NTSPP:NTS)),ETA);
                        
                        %% FORCES FILE
                        
                        fid_force = fopen(strcat(folder_name,'/forces.dat'),'w');
                        fprintf(fid_force,'t y_dot alpha_eff CD CL CM CPP CPR\n');
                        for k = NTS-NTSPP:NTS % PRINT ONLY LAST PERIOD
                            fprintf(fid_force,"%f %f %f %f %f %f %f %f\n",k*dt/(1/freq)-(NPs-1),PLUNG_VEL(k),ALPHA_EFF(k),CD(k),CL(k),-CM(k),-CD(k),CPR(k));
                        end
                        fclose(fid_force);
                                            
                    end
                end
            end
        end
    end
end

%% END
WarnWave = [sin(1:0.25:500), sin(1:0.5:500), sin(1:1.0:500)];
Audio = audioplayer(WarnWave, 11025);
play(Audio);

toc
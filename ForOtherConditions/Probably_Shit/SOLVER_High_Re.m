%  UNSTEADY HSPM by Komasho (Based on TENG)
%  Last Update: June 19th, 2022 by Emanuel Camacho

clc; close all; %clear all;

%% SIMULATION PARAMETERS

% AERODYNAMIC CHORD
chord = 1;
% ANGLE OF ATTACK
AoA = 0*pi/180;
% REYNOLDS NUMBER
Re = [2e5];
% NONDIMENSIONAL AMPLITUDE
NonDimAmpl = [0.0]';
% REDUCED FREQUENCY
RedFreq = [0.01]';
% MEAN ANGLE OF ATTACK
%alpha_0 = (10-20)*pi/180;
%alpha_1 = -(0+20)*pi/180;
alpha_0 = 0*pi/180;
alpha_1 = 0*pi/180;
% PITCHING AMPLITUDE
PitchAmp = -[-10]'*pi/180;
PitchAmp_2 = -[0]'*pi/180;
% PIVOT (CHORD PERCENTAGE)
pivot = 0.3;
% ASYMMETRY OF KINEMATICS
zeta = [0.5]';
% ANIMATION ?
animation = 1;

fparstu = fopen('Results/Results_ParametricStudy.dat','w');
fprintf(fparstu,'VARIABLES=Re,h,k,Alpha,Alpha_2,CD,CL,CM\n');

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
                        freq = RedFreq(iRedFreq)*V_mag/(pi*chord);
                        %freq = RedFreq(iRedFreq)*V_mag/(2*pi*chord);
                        % MOTION PERIOD
                        T = 1/freq;
                        
                        fprintf('======= CURRENT CONDITION =======\n');
                        fprintf('k: %2.3f | h: %2.3f | A: %2.3f | A_2: %2.3f\n',RedFreq(iRedFreq),NonDimAmpl(iNonDimAmpl),PitchAmp(iA_alpha)*180/pi,PitchAmp_2(iA_alpha_2)*180/pi);
                        
                        folder_name = sprintf('Results/Re_%6.0f_h%0.3f_k%0.2f_A%0.2f',Re(iRe),NonDimAmpl(iNonDimAmpl),RedFreq(iRedFreq),PitchAmp(iA_alpha)*180/pi);
                        mkdir(folder_name);
                                                
                        %% NUMERICAL PARAMETERS
                        
                        % NUMBER OF SIMULATED PERIODS
                        NPs = 7;
                        % NUMBER OF TIME STEPS PER PERIOD
                        NTSPP = 100;
                        % TIME STEP CALCULATION (IN SECONDS)
                        dt = (1/freq)/NTSPP;
                        % NUMBER OF TIMESTEPS
                        NTS = NPs*NTSPP;
                        % MAX NUMBER OF ITERATIONS PER TIMESTEP
                        NITER = 50;
                        % MAX NUMBER OF VORTICES
                        NVs = 2*NTSPP;
                        % CONVERGENCE REPORT
                        reportConvergence = 5;
                        
                        %% AIRFOIL COORDINATES
                        
                        % ALERT: (COORDINATES MUST ALWAYS BE CREATED/READ ON THE COUNTERCLOCKWISE DIRECTION)
                        % NUMBER OF BOUNDARY POINTS (ODD NUMBER)
                        N_P = 101;
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
                        %x_wd = x;
                        %y_wd = y;
                                                                        
                        % DETERMINE LEADING EDGE COORDINATES
                        for i = 1:N
                            if x(i) > 0 && x(i+1) < 0
                                iSLE = i;
                            elseif x(i) < 0 && x(i+1) > 0
                                iELE = i+1;
                            end
                        end
                        % UNCOMMENT IF THE WHOLE AIRFOIL IS PITCHING
                        iSLE = -1000;
                        iELE = 1000;
                        
                        for i = 1:N+1
                            if i > iSLE && i < iELE
                                x_wd(i) = x(i)*cos(alpha_0)-y(i)*sin(alpha_0);
                                y_wd(i) = x(i)*sin(alpha_0)+y(i)*cos(alpha_0);
                            else
                                x_wd(i) = x(i)*cos(alpha_1)-y(i)*sin(alpha_1);
                                y_wd(i) = x(i)*sin(alpha_1)+y(i)*cos(alpha_1);
                            end
                        end
                        
                        %% RESULTS FILE
                        
                        % PARAMETRIC STUDY
                        % WAKE VORTICES FILES
                        %fid = fopen('Wake_Results.dat','w');
                        fid = fopen(strcat(folder_name,'/wake.dat'),'w');
                        %fid = fopen(sprintf('Results/Wake_h%0.3f_k%0.2f_A%0.2f.dat',NonDimAmpl(iNonDimAmpl),RedFreq(iRedFreq),PitchAmp(iA_alpha)*180/pi),'w');
                        fprintf(fid,'VARIABLES=X,Y,VX,VY,CV\n');
                        % PRESSURE AND VELOCITY FILES
                        fvec = fopen(strcat(folder_name,'/pressure.dat'),'w');
                        %fvec = fopen(sprintf('Results/PRESS_VEL_h%0.3f_k%0.2f_A%0.2f.dat',NonDimAmpl(iNonDimAmpl),RedFreq(iRedFreq),PitchAmp(iA_alpha)*180/pi),'w');
                        fprintf(fvec,'VARIABLES=X,Y,VX,VY,CPX,CPY,CP\n');
                        
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
                        ALPHA_EFF = zeros([NTS,1]);
                        % EXPANSION RATIO (GEOMETRIC SERIES)
                        RATIO = 1.1;
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
                        CP = zeros([N,1]);
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
                        % SOURCE STRENGTH DISTRIBUTION
                        q = zeros([N,1]);
                        % INITIALIZE VORTICITY
                        TAU = 0;
                        TAUM1 = 0;
                        % INITIALIZE TRAILING-EDGE PANEL PARAMETERS
                        UW_old = V_INF(1);
                        VW_old = V_INF(2);
                        UW = V_INF(1);
                        VW = V_INF(2);
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
                        
                        for k = 1:NTS
                            
                            %% GEOMETRIC PARAMETERS OF THE AIRFOIL
                            
                            for i = 1:N
                                % ASYMMETRIC PLUNGING
                                if k >= floor(k/NTSPP)*NTSPP && k < floor(k/NTSPP)*NTSPP+NTSPP/2
                                    dt = zeta*(1/freq)/(NTSPP/2);
                                    HY = AMP*cos(pi*freq/zeta*(dt*k))-AMP;
                                    V_XY_AIRFOIL(i,:) = [0 -pi*freq*AMP/zeta*sin(pi*freq/zeta*(dt*k))];
                                elseif k >= floor(k/NTSPP)*NTSPP+NTSPP/2 && k < floor(k/NTSPP)*NTSPP+NTSPP
                                    dt = (1-zeta)*(1/freq)/(NTSPP/2);
                                    HY = AMP*cos(pi*freq/(1-zeta)*(dt*k))-AMP;
                                    V_XY_AIRFOIL(i,:) = [0 -pi*freq*AMP/(1-zeta)*sin(pi*freq/(1-zeta)*(dt*k))];
                                end
                                % Y-COORDINATE DEVIATION
                                DHY = HY-HYO;
                                % PLUNGING VELOCITY
                                PLUNG_VEL(k) = -V_XY_AIRFOIL(i,2);
                                
                                % PITCHING
                                %V_AL_NODE_AIRFOIL(i,:) = [0 0];
                                %V_AL_AIRFOIL(i,:) = [0 0];
                                %ROT(k) = 0;
                                
                                % MOVE LEADING EDGE
                                if i > iSLE && i < iELE
                                    %OMEGA = -[0 0 -2*pi*freq*AMP_AL*sin(2*pi*freq*dt*k+pi/2)];
                                    OMEGA = -[0 0 -2*pi*freq*AMP_AL*sin(2*pi*freq*dt*k)]; % angle of attack is symetrical
                                    ROT(k) = OMEGA(3);
                                    % ROTATIONAL VELOCITY AT AIRFOIL NODES
                                    vecpos_node = [x(i) y(i) 0];
                                    vel_node = cross(OMEGA,vecpos_node);
                                    V_AL_NODE_AIRFOIL(i,:) = [vel_node(1) vel_node(2)];
                                    % ROTATIONAL VELOCITY AT THE CONTROL POINTS
                                    vecpos = [0.5*(x(i)+x(i+1)) 0.5*(y(i)+y(i+1)) 0];
                                    vel = cross(OMEGA,vecpos);
                                    V_AL_AIRFOIL(i,:) = [vel(1) vel(2)];
                                else
                                    OMEGA = -[0 0 -2*pi*freq*AMP_AL_2*sin(2*pi*freq*dt*k+pi/2)];
                                    ROT(k) = OMEGA(3);
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
                            ALPHA_1 = AMP_AL*sin(2*pi*freq*dt*k);
                            ALPHA_2 = atan(pivot*sin(ALPHA_1)/(pivot*cos(ALPHA_1)+(1-pivot)));
                            ALPHA_EFF(k) = atan(-PLUNG_VEL(k)/V_mag) - ALPHA_2;
                            ALPHA_EFF(k) = ALPHA_EFF(k)*180/pi;
                            
                            % DEFLECTION IN BODY AXIS (ONLY PITCHING)
                            for i = 1:N+1
                                if i > iSLE && i < iELE
                                    x(i) = x_wd(i)*cos(AMP_AL*cos(2*pi*freq*dt*k)-AMP_AL)-y_wd(i)*sin(AMP_AL*cos(2*pi*freq*dt*k)-AMP_AL);
                                    y(i) = x_wd(i)*sin(AMP_AL*cos(2*pi*freq*dt*k)-AMP_AL)+y_wd(i)*cos(AMP_AL*cos(2*pi*freq*dt*k)-AMP_AL);
                                else
                                    x(i) = x_wd(i)*cos(AMP_AL_2*sin(2*pi*freq*dt*k))-y_wd(i)*sin(AMP_AL_2*sin(2*pi*freq*dt*k));
                                    y(i) = x_wd(i)*sin(AMP_AL_2*sin(2*pi*freq*dt*k))+y_wd(i)*cos(AMP_AL_2*sin(2*pi*freq*dt*k));
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
                            
                            %% ITERATIVE PROCEDURE TO OBTAINED WAKE ELEMENT
                            
                            for iter = 1:NITER
                                
                                %% INFLUENCE COEFFICIENTS
                                
                                % A PANELS INFLUENCE
                                [AN,AT] = calc_AN_AT(x,y,N);
                                % A SHED VORTICITY INFLUENCE
                                [AYW_J,AXW_J] = calc_AXW_J_AYW_J(x,y,N);
                                % A CORE VORTEX INFLUENCE
                                if NM >= 1
                                    [AYH_J,AXH_J] = calc_AXH_J_AYH_J(x,y,XM,YM,N,NM);
                                end
                                % B PANELS INFLUENCE
                                [BN,BT] = calc_BN_BT(x,y,N);
                                % B SHED VORTICITY INFLUENCE
                                [BNI_W,BTI_W] = calc_BNI_W_BTI_W(x,y,N);
                                [BYW_J,BXW_J] = calc_BXW_J_BYW_J(x,y,N);
                                if NM >= 1
                                    [BYH_W,BXH_W] = calc_BXH_W_BYH_W(x,y,XM,YM,N,NM);
                                end
                                % B CORE VORTEX INFLUENCE
                                if NM >= 1
                                    [BYH_J,BXH_J] = calc_BXH_J_BYH_J(x,y,XM,YM,N,NM);
                                end
                                % CALCULATE C COEFFICIENTS IF THERE ARE VORTICES
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
                                                               
                                % CHECK CONVERGENCE (ABSOLUTE)
                                %if abs(VW - VW_old)/abs(VW_old) <= 0.01 && abs(UW - UW_old)/abs(UW_old) <= 0.01
                                if abs(VW - VW_old) <= 0.01 && abs(UW - UW_old) <= 0.01
                                    break;
                                end
                                
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
                            alpha = AoA;
                            
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
                                CX = CX + CP(i)*(y(i+1)-y(i)); %ARE THE SIGNS OK?
                            end
                            
                            %% Y FORCE COEFFICIENTS
                            CY = 0;
                            for i = 1:N
                                CY = CY - CP(i)*(x(i+1)-x(i)); %ARE THE SIGNS OK?
                            end
                            
                            %% MOMENT COEFFICIENTS
                            C_M = 0;
                            for i = 1:N
                                if i > iSLE && i < iELE % CALCULATE THE PITCHING COEFFCIIENT OF MOVABLE PANELS
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
                                           
                            if animation == 1
                                plot(x,y,'black','LineWidth',1);
                                hold on;
                                quiver(CPX,CPY,cprx,cpry,1,'black');
%                                 hold on;
%                                 scatter(XM,YM,10,CV,'filled');
%                                 hold on;
%                                 quiver(0,0,CD(k),0,1,'blue','filled');
%                                 hold on;
%                                 quiver(0,0,0,CL(k),1,'blue','filled');
                                hold off;
                                daspect([1 1 1]);
                                pbaspect([2 1 1]);
                                pause(0.01);
                                %scatter(x,y,10);
                                %plot(XZ,YZ,'o');
                            end
                            
                            %% DATA FOR TECPLOT
                            
                            % PRESURE AND VELOCITY DISTRIBUTION
                            fprintf(fvec,"ZONE T='RESULTS' I=%d\n",N);
                            for airfoil_c = 1:N
                                fprintf(fvec,"%f %f %f %f %f %f %f\n",CPX(airfoil_c),CPY(airfoil_c),vtx(airfoil_c),vty(airfoil_c),cprx(airfoil_c),cpry(airfoil_c),CP(airfoil_c));
                            end
                            
                            % VORTEX
                            fprintf(fid,"ZONE T='RESULTS_h%0.2f_k%0.2f' I=%d\n",NonDimAmpl(iNonDimAmpl),RedFreq(iRedFreq),NM+N+1);
                            for airfoil_c = 1:N+1
                                fprintf(fid,"%f %f %f %f %f\n",x(airfoil_c),y(airfoil_c),V_AL_NODE_AIRFOIL(airfoil_c,1),V_AL_NODE_AIRFOIL(airfoil_c,2),0);
                            end
                            
                            % WAKE
                            for vor = 1:NM
                                fprintf(fid,"%f %f %f %f %f\n",XM(vor),YM(vor),UM(vor),VM(vor),CV(vor));
                            end
                            
                            if rem(k,NTSPP) == 0 && floor(k/NTSPP) > 1
                                %% MEAN DRAG
                                meanCD = mean(CD(k-NTSPP:k));
                                %% MEAN LIFT
                                meanCL = mean(CL(k-NTSPP:k));
                                %% MEAN MOMENT
                                meanCM = mean(CM(k-NTSPP:k));
                                
                                fprintf('Timestep: %3d | meanCD: %2.2f | meanCL: %2.2f | meanCM: %2.2f\n',k,meanCD,meanCL,meanCM);
                                
                            end
                        end
                        
                        fclose(fid);
                        fclose(fvec);
                        
                        %% FORCES DISPLAY
                        
                        fprintf('======= RESULTS =======\n');
                        fprintf('Mean Drag Coefficient (CD) : %2.2f\n',mean(CD(NTS-NTSPP:NTS)));
                        fprintf('Mean Lift Coefficient (CL) : %2.2f\n',mean(CL(NTS-NTSPP:NTS)));
                        fprintf('Mean Moment Coefficient (CM) : %2.2f\n\n',mean(CM(NTS-NTSPP:NTS)));
                        
                        fprintf(fparstu,'%2.5f %2.5f %2.5f %2.5f %2.5f %2.5f %2.5f %2.5f\n',Re(iRe),NonDimAmpl(iNonDimAmpl),RedFreq(iRedFreq),PitchAmp(iA_alpha)*180/pi,PitchAmp_2(iA_alpha_2)*180/pi,mean(CD(NTS-NTSPP:NTS)),mean(CL(NTS-NTSPP:NTS)),mean(CM(NTS-NTSPP:NTS)));
                        
                        %% FORCES FILE
                        
                        fid_force = fopen(strcat(folder_name,'/forces.dat'),'w');
                        for k = NTS-4*NTSPP:NTS % PRINT ONLY LAST PERIOD
                            fprintf(fid_force,"%f %f %f %f %f %f\n",k*dt/(1/freq),PLUNG_VEL(k),ALPHA_EFF(k),CD(k),CL(k),CM(k));
                        end
                        fclose(fid_force);
                        
                        %% PLOTS
                                               
                        figure(1);
%                         plot(x,y,'black','LineWidth',1);
%                         hold on;
%                         pause(0.01);
                        %scatter(x,y,10);
%                         plot(XZ,YZ,'o');
                        %hold on;
                        quiver(CPX,CPY,cprx,cpry,5,'black');
                        %scatter(XM,YM,10,CV,'filled');
                        hold on;
                        quiver(0,0,CD(k),0,0.5,'blue','filled');
                        hold on;
                        quiver(0,0,0,CL(k),0.5,'blue','filled');
                        hold on;
                        quiver(0,0,CX,0,0.5,'red','filled');
                        hold on;
                        quiver(0,0,0,CY,0.5,'red','filled');
                        hold off;
                        daspect([1 1 1]);
                        pbaspect([2 1 1]);

                        figure(2);
                        subplot(3,1,1);
                        plot((1:k)*dt/T,CD(1:k),'black','LineWidth',1);
                        %hold on;
                        %plot((1:k)*dt,CD(1:k)+0.455/log10(Re)^2.58*L/chord,'blue','LineWidth',1);
                        xlabel('$t/T$','interpreter','latex','fontsize',15);
                        ylabel('$C_d$','interpreter','latex','fontsize',15,'rotation',0);

                        subplot(3,1,2);
                        plot((1:k)*dt/T,CL(1:k),'black','LineWidth',1);
                        xlabel('$t/T$','interpreter','latex','fontsize',15);
                        ylabel('$C_l$','interpreter','latex','fontsize',15,'rotation',0);

                        subplot(3,1,3);
                        plot((1:k)*dt/T,CM(1:k),'black','LineWidth',1);
                        xlabel('$t/T$','interpreter','latex','fontsize',15);
                        ylabel('$C_m$','interpreter','latex','fontsize',15,'rotation',0);
                        
                        figure(3);
                        subplot(1,3,1);
                        plot(ALPHA_EFF(NTS-NTSPP:NTS),CD(NTS-NTSPP:NTS),'*');
                        xlabel('$\alpha_{eff}$','interpreter','latex','fontsize',15);
                        ylabel('$C_d$','interpreter','latex','fontsize',15,'rotation',0);

                        subplot(1,3,2);
                        plot(ALPHA_EFF(NTS-NTSPP:NTS),CL(NTS-NTSPP:NTS),'*');
                        xlabel('$\alpha_{eff}$','interpreter','latex','fontsize',15);
                        ylabel('$C_l$','interpreter','latex','fontsize',15,'rotation',0);

                        subplot(1,3,3);
                        plot(ALPHA_EFF(NTS-NTSPP:NTS),CM(NTS-NTSPP:NTS),'*');
                        xlabel('$\alpha_{eff}$','interpreter','latex','fontsize',15);
                        ylabel('$C_m$','interpreter','latex','fontsize',15,'rotation',0);
                        
                    end
                end
            end
        end
    end
end

fclose(fparstu);

%% END
WarnWave = [sin(1:0.25:500), sin(1:0.5:500), sin(1:1.0:500)];
Audio = audioplayer(WarnWave, 11025);
play(Audio);
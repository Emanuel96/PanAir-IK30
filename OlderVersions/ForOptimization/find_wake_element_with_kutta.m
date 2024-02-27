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
        [BYH_W,BXH_W] = calc_BXH_W_BYH_W(x,y,XM,YM,N,NM,NTS);
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
    C3 = D21^2-D2N^2+2*L*TAUM1/dt; %CHECH THIS SHIT IMMEDIATLY
    RADI = sqrt(C2^2-4*C1*C3);
    
    % AIRFOIL VORTICITY
    
    if C1 == 0
        fprintf('Kutta condition may be compromised. Check C1 parameter.\n');
        C1 = 1;
    end
    
    %TAU = roots([C1 C2 C3]); %(-C2-RADI)/(2*C1);
    %TAU = TAU(2);
    
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
    
    [UW,VW] = UW_VW(AXW_J,AYW_J,q,BXW_J,BYW_J,TAU,V_INF,NM,CXW_M,CYW_M,CV);
    
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
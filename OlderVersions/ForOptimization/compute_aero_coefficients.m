% TEN CHORDS AWAY FROM LEADING EDGE
% EXPANSION RATIO (GEOMETRIC SERIES)
RATIO = 1.1;
AX = (10*(1-RATIO))/(1-RATIO^(N_PHI+1));

i_LE = (N_P-1)/2+1;
XZ(1) = x(i_LE);
YZ(1) = y(i_LE);
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
    if onlyLE == 1
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

CD(k) = CD(k)/chord;
CL(k) = CL(k)/chord;
CM(k) = CM(k)/chord;
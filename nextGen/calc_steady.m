
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

% A PANELS INFLUENCE
[AN,AT] = calc_AN_AT(x,y,N);
% B PANELS INFLUENCE
[BN,BT] = calc_BN_BT(x,y,N);

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

fprintf(fvec,"ZONE T='RESULTS' I=%d\n",N);
for airfoil_c = 1:N
    fprintf(fvec,"%f %f %f %f %f %f %f\n",CPX(airfoil_c),CPY(airfoil_c),vtx(airfoil_c),vty(airfoil_c),cprx(airfoil_c),cpry(airfoil_c),CP(airfoil_c));
end

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
for i = i_LE:N+1
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

%% PLOT STEADY SOLUTION

% if NPs == 0
%     plot(x(1:N),y(1:N),'LineWidth',1);
%     hold on;
%     quiver(CPX,CPY,cprx,cpry,1,'black');
%     daspect([1 1 1]);
%     return;
% else
%     %                         plot(x(1:N),y(1:N),'LineWidth',1);
%     %                         hold on;
%     %                         quiver(CPX,CPY,cprx,cpry,1,'black');
%     %                         daspect([1 1 1]);
%     %                         pause(0.01);
% end

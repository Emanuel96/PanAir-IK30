
%% KINEMATICS FILE

for i = 1:N
    
    % VERTICAL MOTION
    HY = y_pos(k);
    V_XY_AIRFOIL(i,:) = [0 y_dot(k)];
    % Y-COORDINATE DEVIATION
    DHY = HY-HYO;
    % PLUNGING VELOCITY
    PLUNG_VEL(k) = -V_XY_AIRFOIL(i,2);
    
    %% TRADITIONAL ROTATION
            
    if i > iSLE && i < iELE
        OMEGA = [0 0 alpha_dot(k)];
        ROT(k) = -OMEGA(3); %IS THE SIGN RIGHT? DANGEROUS
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

% DEFLECTION IN BODY AXIS (ONLY PITCHING)
for i = 1:N+1
    if onlyLE == 0
        pos = alpha_pos(k);
        x(i) = x_wd(i)*cos(pos)-y_wd(i)*sin(pos);
        y(i) = x_wd(i)*sin(pos)+y_wd(i)*cos(pos);
    elseif onlyLE == 1
        pos = AMP_AL/(1+exp(smoother*(x_wd(i)/chord)))*alpha_pos(k);
        x(i) = x_wd(i)*cos(pos)-y_wd(i)*sin(pos);
        y(i) = x_wd(i)*sin(pos)+y_wd(i)*cos(pos);
    end
end

% EFFECTIVE ANGLE OF ATTACK

ALPHA_1(k) = -(alpha_pos(k)-AoA);

% if onlyLE == 1
%     ALPHA_2(k) = AoA; %POSITIVE UP
% elseif onlyLE == 0
%     ALPHA_2(k) = ALPHA_1(k);
% end

%ALPHA_3(k) = atan((pivot*sin(ALPHA_1(k))+(1-pivot)*sin(ALPHA_2(k)))/(pivot*cos(ALPHA_1(k))+(1-pivot)*cos(ALPHA_2(k)))); %ERROR

ALPHA_EFF(k) = atan(-PLUNG_VEL(k)/V_mag) + ALPHA_1(k);
ALPHA_EFF(k) = ALPHA_EFF(k)*180/pi;


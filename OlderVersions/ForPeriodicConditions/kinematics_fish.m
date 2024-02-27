
%% KINEMATICS FILE

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
    
end

%% FISHY DEFORMATION

for i = 1:N+1
    x(i) = x_wd(i);
    y(i) = y_wd(i)+0.1*x_wd(i)^2*sin(2*pi*freq*dt*k);%0.1*x_wd(i)^2*sin(2*pi*x_wd(i)/1.1-2*pi*freq*dt*k);
    x_m1 = x_wd(i);
    x_p1 = x_wd(i);
    y_m1 = y_wd(i)+0.1*x_wd(i)^2*sin(2*pi*freq*(dt*k-dt));%0.1*x_wd(i)^2*sin(2*pi*x_wd(i)/1.1-2*pi*freq*(dt*k-dt));
    y_p1 = y_wd(i)+0.1*x_wd(i)^2*sin(2*pi*freq*(dt*k+dt));%0.1*x_wd(i)^2*sin(2*pi*x_wd(i)/1.1-2*pi*freq*(dt*k+dt));
    vx = (x_p1-x_m1)/(2*dt);
    vy = (y_p1-y_m1)/(2*dt);
    V_AL_NODE_AIRFOIL(i,:) = [vx vy]; %IS THE SIGN RIGHT? DANGEROUS
end

for i = 1:N
    V_AL_AIRFOIL(i,:) = 0.5*(V_AL_NODE_AIRFOIL(i,:)+V_AL_NODE_AIRFOIL(i+1,:));
end

% EFFECTIVE ANGLE OF ATTACK
ALPHA_1(k) = -(AMP_AL*sin(2*pi*freq*dt*k+phase)-AoA); %POSITIVE UP
%ALPHA_1(k) = -(atan(-2*pi*freq*AMP/V_mag*sin(2*pi*freq*dt*k))-AoA);
if onlyLE == 1
    ALPHA_2(k) = AoA; %POSITIVE UP
elseif onlyLE == 0
    ALPHA_2(k) = ALPHA_1(k);
end
ALPHA_3(k) = atan((pivot*sin(ALPHA_1(k))+(1-pivot)*sin(ALPHA_2(k)))/(pivot*cos(ALPHA_1(k))+(1-pivot)*cos(ALPHA_2(k)))); %ERROR
ALPHA_EFF(k) = atan(-PLUNG_VEL(k)/V_mag) + ALPHA_3(k);
ALPHA_EFF(k) = ALPHA_EFF(k)*180/pi;

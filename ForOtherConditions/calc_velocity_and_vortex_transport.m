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
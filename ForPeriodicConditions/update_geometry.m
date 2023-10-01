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

%% PLUNGING
% WITH LINEAR INCREASING EFFECTIVE AOA
%y_dot = tan(pi/180*0.00*[1:NTS]*dt.*(1-exp(-9*[1:NTS]*dt)));
y_dot = tan(pi/180*0*alpha_rate*[1:NTS]*dt*V_mag);
y_pos = cumtrapz([1:NTS]*dt,y_dot);
% WITH SINUSOIDAL EFFECTIVE AOA
%y_dot = tan(pi/180*0*sin(2*pi*1/DURATION*[1:NTS]*dt)*V_mag);
%y_pos = cumtrapz([1:NTS]*dt,y_dot);

%% PITCHING
% WITH LINEAR INCREASING EFFECTIVE AOA
%alpha_dot = -pi/180*0.05*ones(NTS,1).*(1-exp(-9*[1:NTS]*dt));
alpha_dot = -pi/180*1*alpha_rate*ones(NTS,1);
alpha_pos = cumtrapz([1:NTS]*dt,alpha_dot);
% WITH SINUSOIDAL EFFECTIVE AOA
%alpha_pos = -pi/180*10*sin(2*pi*1/DURATION*[1:NTS]*dt);
%alpha_dot = gradient(alpha_pos,dt);


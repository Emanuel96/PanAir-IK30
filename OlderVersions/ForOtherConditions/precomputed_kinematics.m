
if motion == 1
    cv = 1;
elseif motion == 2
    cv = 0;
end

%% PLUNGING
% WITH LINEAR INCREASING EFFECTIVE AOA
%y_dot = tan(pi/180*cv*alpha_rate*[1:NTS]*dt*V_mag);
%y_pos = cumtrapz([1:NTS]*dt,y_dot);

% WITH SINUSOIDAL EFFECTIVE AOA
%y_dot = tan(pi/180*alpha_amp*cv*sin(2*pi*1/DURATION*[1:NTS]*dt)*V_mag);
%y_pos = cumtrapz([1:NTS]*dt,y_dot);

% WITH SIGMOID
y_dot = tan(pi/180*alpha_amp*cv./(1+exp(-0.5*([1:NTS]*dt-DURATION/2))))*V_mag;
y_pos = cumtrapz([1:NTS]*dt,y_dot);

%plot(y_dot)

%% PITCHING
% WITH LINEAR INCREASING EFFECTIVE AOA
%alpha_dot = -pi/180*(1-cv)*alpha_rate*ones(NTS,1);
%alpha_pos = cumtrapz([1:NTS]*dt,alpha_dot);

% WITH SINUSOIDAL EFFECTIVE AOA
%alpha_pos = -pi/180*alpha_amp*(1-cv)*sin(2*pi*1/DURATION*[1:NTS]*dt);
%alpha_dot = gradient(alpha_pos,dt);

% WITH SIGMOID
alpha_pos = -pi/180*alpha_amp*(1-cv)./(1+exp(-0.5*([1:NTS]*dt-DURATION/2)));
alpha_dot = gradient(alpha_pos,dt);


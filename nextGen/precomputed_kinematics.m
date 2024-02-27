
if read_kinematics == 0
    %% PLUNGING
    % SINUSOIDAL PLUNGING
    y_pos = +(AMP*cos(2*pi*1/NTSPP*[1:NTS])-AMP); % IS - CORRECT?
    y_dot = gradient(y_pos,dt);
    y_dot(NTS) = 0;
    
    %% PITCHING
    % SINUSOIDAL PITCHING
    alpha_pos = AMP_AL*sin(2*pi*1/NTSPP*[1:NTS]);
    alpha_dot = gradient(alpha_pos,dt);
else
    
    %% READ KINEMATICS FROM THE kinematics FOLDER
    
    
    
end

%% PLUNGING
% WITH LINEAR INCREASING EFFECTIVE AOA
%y_dot = tan(pi/180*0.00*[1:NTS]*dt.*(1-exp(-9*[1:NTS]*dt)));
%y_pos = cumtrapz([1:NTS]*dt,y_dot);

%% PITCHING
% WITH LINEAR INCREASING EFFECTIVE AOA
%alpha_dot = -pi/180*0.05*ones(NTS,1).*(1-exp(-9*[1:NTS]*dt));
%alpha_pos = cumtrapz([1:NTS]*dt,alpha_dot);

%% RANDOM PITCHING (FOR DATA GENERATION)
% random_alpha_pos = (30*rand(NTS,1)-15)*pi/180;
% f_sample = 10; % [Hz]
% f_c = 1; % 30/(2*pi);
% [bf,af] = butter(10,f_c/(f_sample/2));
% %% FILTERED DATA (BUTTERWORTH FILTER)
% random_alpha_pos_f = filtfilt(bf,af,random_alpha_pos);
% alpha_pos = random_alpha_pos_f.*([1./(1+exp(-0.1.*([1:NTS]-100)))]'+[1./(1+exp(0.1.*([1:NTS]-400)))]'-1);
% alpha_dot = gradient(alpha_pos,dt);
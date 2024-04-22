
%% CONTROL BLOCK

%error(k) = 0.2/(1+exp(-0.05*(k-500))) - CL(k);

%SP = 0.1*sin(2*pi*freq*k*dt);
SP = 0.2/(1+exp(-0.01*(k-250)));

error(k) = SP - CL(k);
int = int + error(k)*dt;

alpha_pos(k+1) = alpha_pos(k) - (0.01*error(k)+0.000*int);
alpha_dot(k+1) = (alpha_pos(k+1)-alpha_pos(k))/dt;

%y_pos(k+1) = y_pos(k) - (0.002*error(k)+0.000*int);
%y_dot(k+1) = (y_pos(k+1)-y_pos(k))/dt;

% if rem(k,5) == 0
% plot(k,SP,'b+');
% hold on;
% plot(k,CL(k),'k+');
% hold on;
% pause(0.01);
% end
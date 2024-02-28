
%% CONTROL BLOCK

%error(k) = 0.2/(1+exp(-0.05*(k-500))) - CL(k);

SP = 0.1*sin(2*pi*0.1*k*dt);

error(k) = SP - CL(k);

int = int + error(k)*dt;

alpha_pos(k+1) = alpha_pos(k) - (0.05*error(k)+0.001*int);

alpha_dot(k+1) = (alpha_pos(k+1)-alpha_pos(k))/dt;

%meanAoA = meanAoA - 0.01*error(k);

plot(k,SP,'b+');
hold on;
plot(k,CL(k),'k+');
hold on;
pause(0.01);
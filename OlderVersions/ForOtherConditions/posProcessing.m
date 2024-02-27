clear clc

%% CONDITIONS

file = 'Results/Re_10000/Re_10000_forces_pitching.dat';
file1 = 'Results/Re_10000/Re_10000_forces_plunging.dat';

result = importdata(file);
result1 = importdata(file1);

time = result.data(:,1);

alpha_eff = result.data(:,3);
drag = result.data(:,4);
lift = result.data(:,5);
moment = result.data(:,6);

alpha_eff1 = result1.data(:,3);

drag1 = result1.data(:,5).*sind(alpha_eff1)+result1.data(:,4).*cosd(alpha_eff1);
lift1 = result1.data(:,5).*cosd(alpha_eff1)-result1.data(:,4).*sind(alpha_eff1);
%drag1 = result1.data(:,4);
%lift1 = result1.data(:,5);
moment1 = result1.data(:,6);

subplot(2,3,1)
plot(time,drag,'b','Linewidth',2);
hold on;
plot(time,drag1,'r');
ylabel('$C_D$','Interpreter','latex','Rotation',0);
xlabel('$t (s)$','Interpreter','latex','Rotation',0);
subplot(2,3,2)
plot(time,lift,'b','Linewidth',2);
hold on;
plot(time,lift1,'r');
ylabel('$C_L$','Interpreter','latex','Rotation',0);
xlabel('$t (s)$','Interpreter','latex','Rotation',0);
subplot(2,3,3)
plot(time,moment,'b','Linewidth',2);
hold on;
plot(time,moment1,'r');
ylabel('$C_M$','Interpreter','latex','Rotation',0);
xlabel('$t (s)$','Interpreter','latex','Rotation',0);
subplot(2,3,4)
plot(alpha_eff,drag,'b','Linewidth',2);
hold on;
plot(alpha_eff1,drag1,'r');
ylabel('$C_D$','Interpreter','latex','Rotation',0);
xlabel('$\alpha [^\circ]$','Interpreter','latex','Rotation',0);
subplot(2,3,5)
plot(alpha_eff,lift,'b','Linewidth',2);
hold on;
plot(alpha_eff1,lift1,'r');
ylabel('$C_L$','Interpreter','latex','Rotation',0);
xlabel('$\alpha [^\circ]$','Interpreter','latex','Rotation',0);
subplot(2,3,6)
plot(alpha_eff,moment,'b','Linewidth',2);
hold on;
plot(alpha_eff1,moment1,'r');
ylabel('$C_M$','Interpreter','latex','Rotation',0);
xlabel('$\alpha [^\circ]$','Interpreter','latex','Rotation',0);

%% END
WarnWave = [sin(1:0.25:500), sin(1:0.5:500), sin(1:1.0:500)];
Audio = audioplayer(WarnWave, 11025);
play(Audio);

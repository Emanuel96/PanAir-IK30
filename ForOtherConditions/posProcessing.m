clear clc

%% CONDITIONS

file = 'Results/Re_10000/Re_10000_forces_pitching.dat';
file1 = 'Results/Re_10000/Re_10000_forces_plunging.dat';

result = importdata(file);
result1 = importdata(file1);

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

subplot(1,3,1)
plot(alpha_eff,drag,'b');
hold on;
plot(alpha_eff1,drag1,'r');
subplot(1,3,2)
plot(alpha_eff,lift,'b');
hold on;
plot(alpha_eff1,lift1,'r');
subplot(1,3,3)
plot(alpha_eff,moment,'b');
hold on;
plot(alpha_eff1,moment1,'r');

%% END
WarnWave = [sin(1:0.25:500), sin(1:0.5:500), sin(1:1.0:500)];
Audio = audioplayer(WarnWave, 11025);
play(Audio);

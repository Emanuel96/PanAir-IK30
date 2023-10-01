%% UNSTEADY HSPM by Komasho (Based on TENG)
%  Last Update: September 19th, 2023 by Emanuel Camacho

%% !!! CAUTION !!!
% PLEASE BE CAREFUL ABOUT THE REDUCED FREQUENCY DEFINITION
% PLEASE VERIFY ALL GOVERNING PARAMETERS

%% CLEAR EVERYTHING FIRST
clc; close all; clear all;

%% SIMULATION CONDITIONS
% AERODYNAMIC CHORD
chord = 1;
% ANGLE OF ATTACK
meanAoA = 0*pi/180;
%meanAoA = [0 2.5 5 7.5 10]'*pi/180;
% REYNOLDS NUMBER
Re = 1e4;
%Re = [1e4]';

% NONDIMENSIONAL AMPLITUDE
%NonDimAmpl = 0.0;
NonDimAmpl = [0.05]';
% REDUCED FREQUENCY
%RedFreq = 1;
RedFreq = [0.25]';
% PITCHING AMPLITUDE
PitchAmp = -0*pi/180;
%PitchAmp = -[0 2.5 5 7.5 10]'*pi/180;
% PHASE
phase = 0*pi/180;
% PIVOT (CHORD PERCENTAGE)
pivot = 0.3;
% ASYMMETRY OF KINEMATICS
zeta = 0.5;

%% START CALCULATIONS

for iRe = 1:size(Re)
    for iAoA = 1:size(meanAoA)
        for iNonDimAmpl = 1:size(NonDimAmpl)
            for iRedFreq = 1:size(RedFreq)
                for iA_alpha = 1:size(PitchAmp)
                    for iZeta = 1:size(zeta)
                        
                        % MEAN ANGLE OF ATTACK
                        condition = sprintf('Re_%0.0f_h%0.3f_k%0.2f_A%0.2f',Re(iRe),NonDimAmpl(iNonDimAmpl),RedFreq(iRedFreq),-PitchAmp(iA_alpha)*180/pi);
                        file = sprintf(strcat('Results/',condition,'/',condition,'_forces.dat'));
                        
                        data = importdata(file);
                                            
                        subplot(2,1,1)
                        plot(data.data(:,1),data.data(:,4));
                        xlabel('$t/T$','interpreter','latex','fontsize',15);
                        ylabel('$C_d$','interpreter','latex','fontsize',15,'rotation',0);
                        hold on;
                        subplot(2,1,2)
                        plot(data.data(:,1),data.data(:,5));
                        xlabel('$t/T$','interpreter','latex','fontsize',15);
                        ylabel('$C_l$','interpreter','latex','fontsize',15,'rotation',0);
                        hold on;
                        
                    end
                end
            end
        end
    end
end

%% END
WarnWave = [sin(1:0.25:500), sin(1:0.5:500), sin(1:1.0:500)];
Audio = audioplayer(WarnWave, 11025);
play(Audio);
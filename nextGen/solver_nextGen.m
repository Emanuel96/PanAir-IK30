%% UNSTEADY HSPM by Komasho (Based on TENG)
%  Last Update: September 26th, 2023 by Emanuel Camacho

%% !!! CAUTION !!!
% PLEASE BE CAREFUL ABOUT THE REDUCED FREQUENCY DEFINITION
% PLEASE VERIFY ALL GOVERNING PARAMETERS

%% CLEAR EVERYTHING FIRST
clc;
close all;
%clear all;

%% CODE CONTROLS
code_parameters();
%% SIMULATION CONDITIONS
conditions();

%% START CALCULATIONS

for iRe = 1:size(Re)
    for iAoA = 1:size(meanAoA)
        for iNonDimAmpl = 1:size(NonDimAmpl)
            for iRedFreq = 1:size(RedFreq)
                for iA_alpha = 1:size(PitchAmp)
                    for iZeta = 1:size(zeta)
                        
                        % MEAN ANGLE OF ATTACK
                        AoA = meanAoA(iAoA);
                        % AIRFOIL VELOCITY AND MAGNITUDE
                        V_INF = [Re(iRe)*vis/(chord*rho)*cos(AoA) Re(iRe)*vis/(chord*rho)*sin(AoA)];
                        V_mag = (V_INF(1)^2+V_INF(2)^2)^0.5;
                        % MOTION AMPLITUDE
                        AMP = NonDimAmpl(iNonDimAmpl)*chord;
                        % PITCHING AMPLITUDE
                        AMP_AL = PitchAmp(iA_alpha);
                        % MOTION FREQUENCY
                        freq = RedFreq(iRedFreq)*V_mag/(2*pi*chord);
                        % MOTION PERIOD
                        T = 1/freq;
                        % TIME STEP CALCULATION (IN SECONDS)
                        dt = (1/freq)/NTSPP;
                        
                        fprintf('======= CURRENT CONDITION =======\n');
                        fprintf('Re: %0.2f (U = %f m/s) | k: %2.3f (f = %f Hz) | h: %2.3f | A: %2.3f\n',Re(iRe),V_mag,RedFreq(iRedFreq),freq,NonDimAmpl(iNonDimAmpl),-PitchAmp(iA_alpha)*180/pi);
                        folder_name = sprintf('Results/Re_%0.0f_h%0.3f_k%0.2f_A%0.2f',Re(iRe),NonDimAmpl(iNonDimAmpl),RedFreq(iRedFreq),-PitchAmp(iA_alpha)*180/pi);
                        
                        mkdir(folder_name);
                                                
                        %% AIRFOIL DATA (COORDINATES)
                        geometry();
                        %% RESULTS FILE
                        prepare_results_file();
                        %% INITIALIZE VARIABLES
                        init();
                        %% FIND STEADY SOLUTION
                        calc_steady();
                        %% PRECOMPUTATION OF KINEMATICS
                        precomputed_kinematics();
                        
                        %% START OF TRANSIENT CALCULATION
                        for k = 1:NTS
                            %% KINEMATICS
                            kinematics();
                            %% GEOMETRIC PARAMETERS OF THE AIRFOIL
                            update_geometry();
                            %% STREAM VELOCITY AND VORTEX CONVECTION
                            calc_velocity_and_vortex_transport();
                            %% ITERATIVE PROCEDURE TO CALCULATE WAKE ELEMENT
                            find_wake_element_with_kutta();
                            %% COMPUTE AERODYNAMIC COEFFICIENTS COEFFICIENTS
                            compute_aero_coefficients();
                            %% DEATTACH TRAILING-EDGE PANEL AND CREATE CORE VORTEX
                            shed_vortex();
                            %% DATA FOR TECPLOT
                            data_for_tecplot();
                            %% ANIMATIONS
                            animations();
                            %% SHOW CURRENT STATE
                            simulation_state();
                        end
                        
                        %% FORCES DISPLAY AND FILES
                        display_forces_and_write();
                        %% PLOTS
                        plots();
                        
                        %% CLOSE FILES
                        fclose(fid);
                        fclose(fvec);
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
%% UNSTEADY HSPM by Komasho (Based on TENG)
%  Last Update: September 26th, 2023 by Emanuel Camacho

%% !!! CAUTION !!!
% PLEASE BE CAREFUL ABOUT THE REDUCED FREQUENCY DEFINITION
% PLEASE VERIFY ALL GOVERNING PARAMETERS

%% CLEAR EVERYTHING FIRST
clc;
close all;
clear all;

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
                        
                        %% DIMENSIONAL PARAMETERS
                        simulation();
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
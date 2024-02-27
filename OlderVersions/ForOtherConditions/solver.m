%% UNSTEADY HSPM by Komasho (Based on TENG)
%  Last Update: September 26th, 2023 by Emanuel Camacho

%% !!! CAUTION !!!
% PLEASE BE CAREFUL ABOUT THE REDUCED FREQUENCY DEFINITION
% PLEASE VERIFY ALL GOVERNING PARAMETERS

%% CLEAR EVERYTHING FIRST
clc; close all; clear;

%% SIMULATION CONDITIONS
conditions();

%% START CALCULATIONS
%fparstu = fopen('Results/Results_ParametricStudy.dat','w');
%fprintf(fparstu,'VARIABLES=Re,h,k,Alpha,Alpha_2,CPP,CPR,ETA\n');

for iRe = 1:size(Re)
    for iAoA = 1:size(meanAoA)
        
        % MEAN ANGLE OF ATTACK
        AoA = meanAoA(iAoA);
        % AIRFOIL VELOCITY AND MAGNITUDE
        V_INF = [Re(iRe)*vis/(chord*rho)*cos(AoA) Re(iRe)*vis/(chord*rho)*sin(AoA)];
        V_mag = (V_INF(1)^2+V_INF(2)^2)^0.5;
        
        fprintf('======= CURRENT CONDITION =======\n');
        fprintf('Re: %0.0f (U = %f m/s)\n',Re(iRe),V_mag);
        folder_name = sprintf('Results/Re_%0.0f',Re(iRe));
        
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
            if save_data == 1
                data_for_tecplot();
            end
            %% PROPULSIVE AND REQUIRED POWER COEFFICIENTS
            update_propulsive_coefficients();
            %% PLOT
            quiver(CPX,CPY,cprx,cpry);
            daspect([1 1 1])
            pause(0.01);
        end
        
        fclose(fid);
        fclose(fvec);
        
        %% PROPULSIVE PERFORMANCE
        propulsive_coefficients();
        %% FORCES DISPLAY AND FILES
        display_forces_and_write();
        %% PLOTS
        plots();
        %% NEURAL NETS SECTIONS - DATA GENERATION
        %data_generation_NNs();
    end
end

%fclose(fparstu);

%% END
WarnWave = [sin(1:0.25:500), sin(1:0.5:500), sin(1:1.0:500)];
Audio = audioplayer(WarnWave, 11025);
play(Audio);
% CONDITION NUMBER
%folder_name = sprintf('G:/My Drive/[AEROG] NeuroLift/Results/DataForValidation_plunging/%d',cond_no);
%folder_name = sprintf('Results/dataForNeuralLift',cond_no);
folder_name = sprintf('G:/My Drive/[AEROG] NeuroLift/Data/DataForValidation_pitching/%d',cond_no);
mkdir(folder_name);
% REGIME FILE
regime_file = fopen(strcat(folder_name,'/','regime.dat'),'w');
fprintf(regime_file,"Re h k a mean_AoA\n");
fprintf(regime_file,"%f %f %f %f %f",Re(iRe),NonDimAmpl(iNonDimAmpl),RedFreq(iRedFreq),-PitchAmp(iA_alpha)*180/pi,AoA*180/pi);
fclose(regime_file);
% FORCES FILE
fid_force = fopen(strcat(folder_name,'/','forces.dat'),'w');
fprintf(fid_force,'t y_dot alpha_eff CD CL CM\n');
for k = NTS-2*NTSPP:NTS % PRINT ONLY LAST PERIOD
    fprintf(fid_force,"%f %f %f %f %f %f\n",k*dt/(1/freq)-(NPs-1),PLUNG_VEL(k),ALPHA_EFF(k),CD(k),CL(k),-CM(k));
end
fclose(fid_force);

cond_no = cond_no + 1;
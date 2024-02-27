
%% FORCES FILE

fid_force = fopen(strcat(folder_name,'/',condition,'_forces.dat'),'w');
fprintf(fid_force,'t y_dot alpha_eff CD CL CM\n');
for k = NTS-NTSPP:NTS % PRINT ONLY LAST PERIOD
    fprintf(fid_force,"%f %f %f %f %f %f\n",k*dt/(1/freq),PLUNG_VEL(k),ALPHA_EFF(k),CD(k),CL(k),-CM(k));
end
fclose(fid_force);
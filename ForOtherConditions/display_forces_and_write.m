

% fprintf('======= RESULTS =======\n');
% fprintf('Mean Propulsive Power (CPP) : %2.5f (%2.5f N)\n',mean(CT(NTS-NTSPP:NTS)),1/2*rho*V_mag^2*chord*1*mean(CT(NTS-NTSPP:NTS)));
% fprintf('Mean Power Coefficient (CPR) : %2.5f\n',mean(CPR(NTS-NTSPP:NTS)));
% fprintf('Propulsive Efficiency (ETA) : %2.5f\n\n',ETA);

%fprintf(fparstu,'%2.5f %2.5f %2.5f %2.5f %2.5f %2.5f %2.5f\n',Re(iRe),NonDimAmpl(iNonDimAmpl),RedFreq(iRedFreq),PitchAmp(iA_alpha)*180/pi,mean(CT(NTS-NTSPP:NTS)),mean(CPR(NTS-NTSPP:NTS)),ETA);

%% FORCES FILE

fid_force = fopen(strcat(folder_name,'/',condition,'_forces.dat'),'w');
%fprintf(fid_force,'t y_dot alpha_eff CD CL CM CPP CPR\n');
fprintf(fid_force,'t y_dot alpha_eff CD CL CM\n');
%for k = NTS-NTSPP:NTS % PRINT ONLY LAST PERIOD
for k = 1:NTS % PRINT ONLY LAST PERIOD
    %fprintf(fid_force,"%f %f %f %f %f %f %f %f\n",k*dt/(1/freq)-(NPs-1),PLUNG_VEL(k),ALPHA_EFF(k),CD(k),CL(k),-CM(k),-CD(k),CPR(k));
    fprintf(fid_force,"%f %f %f %f %f %f\n",k*dt,PLUNG_VEL(k),ALPHA_EFF(k),CD(k),CL(k),-CM(k));
end
fclose(fid_force);
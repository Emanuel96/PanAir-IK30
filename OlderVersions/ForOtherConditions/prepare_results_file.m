% PARAMETRIC STUDY
% WAKE VORTICES FILES
%fid = fopen('Wake_Results.dat','w');

condition = sprintf('Re_%0.0f',Re(iRe));

fid = fopen(strcat(folder_name,'/wake.dat'),'w');
%fid = fopen(sprintf('Results/Wake_h%0.3f_k%0.2f_A%0.2f.dat',NonDimAmpl(iNonDimAmpl),RedFreq(iRedFreq),PitchAmp(iA_alpha)*180/pi),'w');
fprintf(fid,'VARIABLES=X,Y,VX,VY,CV\n');
% PRESSURE AND VELOCITY FILES
fvec = fopen(strcat(folder_name,'/pressure.dat'),'w');
%fvec = fopen(sprintf('Results/PRESS_VEL_h%0.3f_k%0.2f_A%0.2f.dat',NonDimAmpl(iNonDimAmpl),RedFreq(iRedFreq),PitchAmp(iA_alpha)*180/pi),'w');
fprintf(fvec,'VARIABLES=X,Y,VX,VY,CPX,CPY,CP\n');
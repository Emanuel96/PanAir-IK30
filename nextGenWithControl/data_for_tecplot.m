
if save_data == 1

vtx = V_t.*cos(theta_t);
vty = V_t.*sin(theta_t);
cprx = -CP.*cos(theta_n);
cpry = -CP.*sin(theta_n);
                        
% PRESURE AND VELOCITY DISTRIBUTION
fprintf(fvec,"ZONE T='RESULTS' I=%d\n",N);
for airfoil_c = 1:N
    fprintf(fvec,"%f %f %f %f %f %f %f\n",CPX(airfoil_c),CPY(airfoil_c),vtx(airfoil_c),vty(airfoil_c),cprx(airfoil_c),cpry(airfoil_c),CP(airfoil_c));
end

% if k == NTS-NTSPP/4
%     press_desc = fopen(strcat(folder_name,'/',condition,'_pressure_075.dat'),'w');
%     for airfoil_c = 1:N
%         fprintf(press_desc,"%f %f %f\n",CPX(airfoil_c),CPY(airfoil_c),CP(airfoil_c));
%     end
%     fclose(press_desc);
% elseif k == NTS
%     press_desc = fopen(strcat(folder_name,'/',condition,'_pressure_100.dat'),'w');
%     for airfoil_c = 1:N
%         fprintf(press_desc,"%f %f %f\n",CPX(airfoil_c),CPY(airfoil_c),CP(airfoil_c));
%     end
%     fclose(press_desc);
% end

% VORTEX
fprintf(fid,"ZONE T='RESULTS_h%0.2f_k%0.2f' I=%d\n",NonDimAmpl(iNonDimAmpl),RedFreq(iRedFreq),NM+N+1);
for airfoil_c = 1:N+1
    fprintf(fid,"%f %f %f %f %f\n",x(airfoil_c),y(airfoil_c),V_AL_NODE_AIRFOIL(airfoil_c,1),V_AL_NODE_AIRFOIL(airfoil_c,2),0);
end

% WAKE
for vor = 1:NM
    fprintf(fid,"%f %f %f %f %f\n",XM(vor),YM(vor),UM(vor),VM(vor),CV(vor));
end

end
if rem(k,NTSPP) == 0 && floor(k/NTSPP) > 1
    CT = -CD;
    %POWER COEFFICIENT
    plungingP = -PLUNG_VEL.*CL/V_mag;
    %CM = -CM; %VERIFY THIS
    pitchingP = -ROT.*(-CM)*chord/V_mag;
    CPR = plungingP+pitchingP;
    
    % MEAN DRAG
    meanCPP = -mean(CD(k-NTSPP:k));
    % MEAN LIFT
    meanCPR = mean(CPR(k-NTSPP:k));
    % EFFICIENCY
    ETA = meanCPP/meanCPR;
    
    fprintf('Timestep: %3d | meanCPP: %2.5f | meanCPR: %2.5f | ETA: %2.5f\n',k,meanCPP,meanCPR,ETA);
    %fprintf('Timestep %3d completed\n',k);
    
end
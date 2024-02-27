

%% OPTIMIZATION FILE

if rem(k,NTSPP) == 0 && floor(k/NTSPP) > 1
    
    DER_1 = (ETA-old_CP)/(AMP_AL-AMP_AL_0);
    
    AMP_AL_0 = AMP_AL;
    
    AMP_AL = AMP_AL_0 + 0.01*DER_1;
    
    disp(AMP_AL*180/pi);

    %DER_0 = DER_1;
    
    old_CP = ETA;
    
    figure(7)
    plot(k,AMP_AL,'*');
    hold on;
    pause(0.01);
    
end

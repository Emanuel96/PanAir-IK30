
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

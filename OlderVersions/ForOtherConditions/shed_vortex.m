
%% SHED VORTEX AND PREPARE FOR NEXT TIME STEP

%% VORTEX CORE VELOCITIES
if NM >= 1
    [UM,VM] = UM_VM(AXH_J,AYH_J,q,BXH_J,BYH_J,BXH_W,BYH_W,TAU,TAUM1,Delta,L,V_INF,NM,CXH_M,CYH_M,CV,k);
end

%% DEATTACH TRAILING-EDGE PANEL AND CREATE CORE VORTEX

NM = NM + 1;

if NM == NVs
    NM = NVs-1;
    UM = circshift(UM,-1);
    VM = circshift(VM,-1);
    XM = circshift(XM,-1);
    YM = circshift(YM,-1);
    CV = circshift(CV,-1);
end

UM(NM) = UW;
VM(NM) = VW;
XM(NM) = x(N+1) + 0.5*Delta*cos(Theta);
YM(NM) = y(N+1) + 0.5*Delta*sin(Theta);
% CORE VORTEX CIRCULATION
CV(NM) = L*(TAUM1-TAU);

%% PREPARATIONS FOR THE NEXT TIMESTEP

HYO = HY;
TAUM1 = TAU;
PHIM1 = PHI;
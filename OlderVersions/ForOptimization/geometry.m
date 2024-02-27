%% AIRFOIL COORDINATES

% ALERT: (COORDINATES MUST ALWAYS BE CREATED/READ ON THE COUNTERCLOCKWISE DIRECTION)
% NUMBER OF BOUNDARY POINTS (ODD NUMBER)
if rem(N_P,2) == 0 % FORCE ODD NUMBER
    N_P = N_P + 1;
end
% NUMBER OF PANELS
N = N_P-1;

%% NACA0012 COORDINATE GENERATOR (WITH COSSINE SPACING)
beta = linspace(0,2*pi,N_P)';
x1 = 0.5*(1-cos(beta+pi));
y1 = 0.594689181*(0.298222773*sqrt(x1)-0.127125232.*x1-0.357907906*x1.^2+0.291984971*x1.^3-0.105174606*x1.^4);
x = (x1-pivot)*chord;
y=([y1(1:(N_P-1)/2);-y1((N_P-1)/2+1:N_P)])*chord;
y=flip(y);
x=round(x,7);
y=round(y,7);

%ALERT: (COORDINATES MUST ALWAYS BE CREATED/READ ON THE COUNTERCLOCKWISE DIRECTION)
%                     SD7003 = readmatrix('AIRFOILS/SD7003_REFINED.dat');
%                     x = (SD7003(:,1)-pivot)*chord;
%                     y = (SD7003(:,2))*chord;
%                     x = flip(x);
%                     y = flip(y);
%                     N_P = size(x,1);
%                     N = N_P-1;

% AIFOIL COORDINATES WITHOUT DEFORMATION (+ROTATION)
x_wd = x;
y_wd = y;

% DETERMINE LEADING EDGE COORDINATES
if onlyLE == 1
    for i = 1:N
        if x(i) > 0 && x(i+1) < 0
            iSLE = i;
        elseif x(i) < 0 && x(i+1) > 0
            iELE = i+1;
        end
    end
else
    iSLE = -1000;
    iELE = 1000;
end
                    
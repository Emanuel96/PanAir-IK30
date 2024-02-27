function [CTF_M] = MATRICES_FOR_FORCES_C(N_PHI,XZ,YZ,alpha,XM,YM,NM)

    %% CTI_M_F
    
    for f = 1:N_PHI
        for m = 1:NM
            
            x_m = 0.5*(XZ(f)+XZ(f+1));
            y_m = 0.5*(YZ(f)+YZ(f+1));
            
            theta = atan2(y_m-YM(m),x_m-XM(m));
            
            % SIN(X-Y) = SIN(X)COS(Y)-COS(X)SIN(Y);
            SINTERM = -sin(alpha)*cos(theta)-(-cos(alpha)*sin(theta));
            
            r = ((x_m-XM(m))^2+(y_m-YM(m))^2)^0.5;

            CTF_M(f,m) = -SINTERM/(2*pi*r);
        end
    end

end


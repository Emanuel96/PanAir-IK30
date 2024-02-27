function [AT_F,BT_F,BTF_W] = MATRICES_FOR_FORCES(N_PHI,N,XZ,YZ,x,y,alpha)

    %% AT_F

    for f = 1:N_PHI
        
        x_m = 0.5*(XZ(f)+XZ(f+1));
        y_m = 0.5*(YZ(f)+YZ(f+1));
        theta = atan2(y(N+2)-y(N+1),x(N+2)-x(N+1));
        % SIN(X-Y) = SIN(X)COS(Y)-COS(X)SIN(Y);
        SINTERM = -sin(alpha)*cos(theta)-(-cos(alpha)*sin(theta));
        % COS(X-Y) = COS(X)COS(Y)+SIN(X)SIN(Y);
        COSTERM = -cos(alpha)*cos(theta)+(-sin(alpha)*sin(theta));
        beta = atan2((y_m-y(N+2))*(x_m-x(N+1))-(x_m-x(N+2))*(y_m-y(N+1)),(x_m-x(N+2))*(x_m-x(N+1))+(y_m-y(N+2))*(y_m-y(N+1)));    
        r = ((x_m-x(N+1))^2+(y_m-y(N+1))^2)^0.5;
        rP1 = ((x_m-x(N+2))^2+(y_m-y(N+2))^2)^0.5;
        BTF_W(f) = 1/(2*pi)*(SINTERM*log(rP1/r)+COSTERM*beta);
        
        for j = 1:N
            %x_m = 0.5*(XZ(f)+XZ(f+1));
            %y_m = 0.5*(YZ(f)+YZ(f+1));

            theta = atan2(y(j+1)-y(j),x(j+1)-x(j));

            % SIN(X-Y) = SIN(X)COS(Y)-COS(X)SIN(Y);
            SINTERM = -sin(alpha)*cos(theta)-(-cos(alpha)*sin(theta));
            % COS(X-Y) = COS(X)COS(Y)+SIN(X)SIN(Y);
            COSTERM = -cos(alpha)*cos(theta)+(-sin(alpha)*sin(theta));

            beta = atan2((y_m-y(j+1))*(x_m-x(j))-(x_m-x(j+1))*(y_m-y(j)),(x_m-x(j+1))*(x_m-x(j))+(y_m-y(j+1))*(y_m-y(j)));   

            r = ((x_m-x(j))^2+(y_m-y(j))^2)^0.5;
            rP1 = ((x_m-x(j+1))^2+(y_m-y(j+1))^2)^0.5;

            AT_F(f,j) = 1/(2*pi)*(SINTERM*beta-COSTERM*log(rP1/r));
            BT_F(f,j) = 1/(2*pi)*(SINTERM*log(rP1/r)+COSTERM*beta);            
        end
                
    end
end


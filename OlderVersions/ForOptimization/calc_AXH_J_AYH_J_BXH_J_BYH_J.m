function [AXH_J,AYH_J,BXH_J,BYH_J] = calc_AXH_J_AYH_J_BXH_J_BYH_J(x,y,XM,YM,N,NM,NTS)
    
    AXH_J = zeros([NTS,N]);
    AYH_J = zeros([NTS,N]);
    BXH_J = zeros([NTS,N]);
    BYH_J = zeros([NTS,N]);

    for h = 1:NM
        
        x_m = XM(h);
        y_m = YM(h);
        theta1 = 0;
        
        for j = 1:N
%             if h == j
%                 AXH_J(h,j) = 0.5;
%                 AYH_J(h,j) = 0;
%             else
                
                theta2 = atan2(y(j+1)-y(j),x(j+1)-x(j));
                beta = atan2((y_m-y(j+1))*(x_m-x(j))-(x_m-x(j+1))*(y_m-y(j)),(x_m-x(j+1))*(x_m-x(j))+(y_m-y(j+1))*(y_m-y(j)));            
                r = ((x_m-x(j))^2+(y_m-y(j))^2)^0.5;
                rP1 = ((x_m-x(j+1))^2+(y_m-y(j+1))^2)^0.5;
                AXH_J(h,j) = 1/(2*pi)*(sin(theta1-theta2)*log(rP1/r)+cos(theta1-theta2)*beta);
                AYH_J(h,j) = 1/(2*pi)*(sin(theta1-theta2)*beta-cos(theta1-theta2)*log(rP1/r));
				
				BXH_J(h,j) = -AYH_J(h,j);
                BYH_J(h,j) = AXH_J(h,j);
%             end
        end
    end
end


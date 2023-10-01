function [BXH_J,BYH_J] = calc_BXH_J_BYH_J(x,y,XM,YM,N,NM)
    for h = 1:NM
        for j = 1:N
%             if h == j
%                 BXH_J(h,j) = 0;
%                 BYH_J(h,j) = 0.5;
%             else
                x_m = XM(h);
                y_m = YM(h);
                theta1 = 0;
                theta2 = atan2(y(j+1)-y(j),x(j+1)-x(j));
                beta = atan2((y_m-y(j+1))*(x_m-x(j))-(x_m-x(j+1))*(y_m-y(j)),(x_m-x(j+1))*(x_m-x(j))+(y_m-y(j+1))*(y_m-y(j)));            
                r = ((x_m-x(j))^2+(y_m-y(j))^2)^0.5;
                rP1 = ((x_m-x(j+1))^2+(y_m-y(j+1))^2)^0.5;
                BXH_J(h,j) = 1/(2*pi)*(cos(theta1-theta2)*log(rP1/r)-sin(theta1-theta2)*beta);
                BYH_J(h,j) = 1/(2*pi)*(cos(theta1-theta2)*beta+sin(theta1-theta2)*log(rP1/r));
%             end
        end
    end
end


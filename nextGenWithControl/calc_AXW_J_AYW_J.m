function [AXW_J,AYW_J] = calc_AXW_J_AYW_J(x,y,N)
    for j = 1:N
        x_m = 0.5*(x(N+1)+x(N+2));
        y_m = 0.5*(y(N+1)+y(N+2));
        theta1 = 0;
        theta2 = atan2(y(j+1)-y(j),x(j+1)-x(j));
        beta = atan2((y_m-y(j+1))*(x_m-x(j))-(x_m-x(j+1))*(y_m-y(j)),(x_m-x(j+1))*(x_m-x(j))+(y_m-y(j+1))*(y_m-y(j))); 
        r = ((x_m-x(j))^2+(y_m-y(j))^2)^0.5;
        rP1 = ((x_m-x(j+1))^2+(y_m-y(j+1))^2)^0.5;
        AXW_J(j) = 1/(2*pi)*(sin(theta1-theta2)*log(rP1/r)+cos(theta1-theta2)*beta);
        AYW_J(j) = 1/(2*pi)*(sin(theta1-theta2)*beta-cos(theta1-theta2)*log(rP1/r));
    end
end


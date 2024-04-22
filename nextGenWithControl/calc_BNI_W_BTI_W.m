function [BNI_W,BTI_W] = calc_BNI_W_BTI_W(x,y,N)
    for i = 1:N
        x_m = 0.5*(x(i)+x(i+1));
        y_m = 0.5*(y(i)+y(i+1));
        theta1 = atan2(y(i+1)-y(i),x(i+1)-x(i));
        theta2 = atan2(y(N+2)-y(N+1),x(N+2)-x(N+1));
        beta = atan2((y_m-y(N+2))*(x_m-x(N+1))-(x_m-x(N+2))*(y_m-y(N+1)),(x_m-x(N+2))*(x_m-x(N+1))+(y_m-y(N+2))*(y_m-y(N+1)));            
        r = ((x_m-x(N+1))^2+(y_m-y(N+1))^2)^0.5;
        rP1 = ((x_m-x(N+2))^2+(y_m-y(N+2))^2)^0.5;
        BNI_W(i) = 1/(2*pi)*(cos(theta1-theta2)*log(rP1/r)-sin(theta1-theta2)*beta);
        BTI_W(i) = 1/(2*pi)*(cos(theta1-theta2)*beta+sin(theta1-theta2)*log(rP1/r));
    end
end


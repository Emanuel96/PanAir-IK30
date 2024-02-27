function [BXH_W,BYH_W] = calc_BXH_W_BYH_W(x,y,XM,YM,N,NM,NTS)

    BXH_W = zeros([NTS,1]);
    BYH_W = zeros([NTS,1]);

    for h = 1:NM
        x_m = XM(h);
        y_m = YM(h);
        theta1 = 0;
        theta2 = atan2(y(N+2)-y(N+1),x(N+2)-x(N+1));
        beta = atan2((y_m-y(N+2))*(x_m-x(N+1))-(x_m-x(N+2))*(y_m-y(N+1)),(x_m-x(N+2))*(x_m-x(N+1))+(y_m-y(N+2))*(y_m-y(N+1)));            
        r = ((x_m-x(N+1))^2+(y_m-y(N+1))^2)^0.5;
        rP1 = ((x_m-x(N+2))^2+(y_m-y(N+2))^2)^0.5;
        BXH_W(h) = 1/(2*pi)*(cos(theta1-theta2)*log(rP1/r)-sin(theta1-theta2)*beta);
        BYH_W(h) = 1/(2*pi)*(cos(theta1-theta2)*beta+sin(theta1-theta2)*log(rP1/r));
    end
end


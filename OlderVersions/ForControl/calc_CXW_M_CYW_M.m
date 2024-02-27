function [CXW_M,CYW_M] = calc_CXW_M_CYW_M(x,y,XM,YM,N,NM)
    for m = 1:NM
        x_m = 0.5*(x(N+1)+x(N+2));
        y_m = 0.5*(y(N+1)+y(N+2));
        theta1 = 0;
        theta2 = atan2(y_m-YM(m),x_m-XM(m));
        r = ((x_m-XM(m))^2+(y_m-YM(m))^2)^0.5;
        CXW_M(m) = -cos(theta1-theta2)/(2*pi*r);
        CYW_M(m) = -sin(theta1-theta2)/(2*pi*r);
    end
end


function [CNI_M,CTI_M] = calc_CNI_M_CTI_M(x,y,XM,YM,N,NM)
    for i = 1:N
        for m = 1:NM
            x_m = 0.5*(x(i)+x(i+1));
            y_m = 0.5*(y(i)+y(i+1));
            theta1 = atan2(y(i+1)-y(i),x(i+1)-x(i));
            theta2 = atan2(y_m-YM(m),x_m-XM(m));
            r = ((x_m-XM(m))^2+(y_m-YM(m))^2)^0.5;
            CNI_M(i,m) = -cos(theta1-theta2)/(2*pi*r);
            CTI_M(i,m) = -sin(theta1-theta2)/(2*pi*r);
        end
    end
end


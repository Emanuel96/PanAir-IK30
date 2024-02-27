function [CXH_M,CYH_M] = calc_CXH_M_CYH_M(XM,YM,NM)
    for h = 1:NM
        for m = 1:NM
            x_m = XM(h);
            y_m = YM(h);
            theta1 = 0;
            theta2 = atan2(y_m-YM(m),x_m-XM(m));
            r = sqrt((x_m-XM(m))^2+(y_m-YM(m))^2);
            if r == 0
                r = 1;
            end
            CXH_M(h,m) = -cos(theta1-theta2)/(2*pi*r);
            CYH_M(h,m) = -sin(theta1-theta2)/(2*pi*r);
        end
    end
end


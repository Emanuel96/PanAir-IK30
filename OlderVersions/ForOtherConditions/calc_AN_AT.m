function [AN,AT] = calc_AN_AT(x,y,N)

    for i = 1:N
        for j = 1:N
            if i == j
                AN(i,j) = 0.5;
                AT(i,j) = 0;
            else
                x_m = 0.5*(x(i)+x(i+1));
                y_m = 0.5*(y(i)+y(i+1));
                theta1 = atan2(y(i+1)-y(i),x(i+1)-x(i));
                theta2 = atan2(y(j+1)-y(j),x(j+1)-x(j));
                beta = atan2((y_m-y(j+1))*(x_m-x(j))-(x_m-x(j+1))*(y_m-y(j)),(x_m-x(j+1))*(x_m-x(j))+(y_m-y(j+1))*(y_m-y(j)));            
                r = ((x_m-x(j))^2+(y_m-y(j))^2)^0.5;
                rP1 = ((x_m-x(j+1))^2+(y_m-y(j+1))^2)^0.5;
                AN(i,j) = 1/(2*pi)*(sin(theta1-theta2)*log(rP1/r)+cos(theta1-theta2)*beta);
                AT(i,j) = 1/(2*pi)*(sin(theta1-theta2)*beta-cos(theta1-theta2)*log(rP1/r));
            end
        end
    end
end
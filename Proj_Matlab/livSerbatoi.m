function x_dot = livSerbatoi(t,x,A,a,k,gamma,g,v)

x_dot = zeros(4,1);

x_dot(1) = (-a(1)*sqrt(2*g*x(1))+a(3)*sqrt(2*g*x(3))+gamma(1)*k(1)*v(1))/A(1);
x_dot(2) = (-a(2)*sqrt(2*g*x(2))+a(4)*sqrt(2*g*x(4))+gamma(2)*k(2)*v(2))/A(2);
x_dot(3) = (-a(3)*sqrt(2*g*x(3)) + (1-gamma(2))*k(2)*v(2))/A(3);
x_dot(4) = (-a(4)*sqrt(2*g*x(4)) + (1-gamma(1))*k(1)*v(1))/A(4);

end
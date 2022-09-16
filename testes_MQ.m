y = [12.2 11.8 11.6 11.6 11.8 12.2 13 14.4 16.2 15.8]'
u = [2.5 2.5 2.5 2.5 2.5 2.23 2.2 2.2 2.21 2.2]'

N_d = length(y);

phi1 = [y(3:N_d-1) y(2:N_d-2) u(1:N_d-3) u(3:N_d-1) u(2:N_d-2)];
theta1 = pinv(phi1)*y(4:N_d)

theta1l=[-5.6841 9.3925 0.7452 -4.2141 -8.3712]';

r1 = phi1*theta1;
Y2 = phi1*theta1l

comp = [y(4:N_d) r1 r2];

MSE1 = mean((y(4:N_d)- r1).^2);
MSE2 = mean((y(4:N_d)- r2).^2);

y2 = [y(1:3);Y2];
 
phi2 = [y2(3:N_d-1) y2(2:N_d-2) u(1:N_d-3) u(3:N_d-1) u(2:N_d-2)];
theta2 = pinv(phi2)*y2(4:N_d)

phi3 = [Y2(3:N_d-3) Y2(2:N_d-3-2) u(4:N_d-3) u(6:N_d-1) u(5:N_d-2)];
theta2 = pinv(phi2)*y2(4:N_d)
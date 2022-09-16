x = 0:199;
u = sin(x)';

% figure
% plot(x,u)
% hold all

u1 = u(1:end/2);
y1 = [u(end/4:end);u(1:end/4-1)];

% plot(u1)
% hold all
% plot(y1)

N = floor(length(u)/2);
ruy_norm = zeros(N,1);
for k=0:N-1
    ruy_norm_T = (u(1:N)'-mean(u))*(y1(1+k:N+k)-mean(y1));
    ruy_norm(k+1) = ruy_norm_T/(2*N+1);
end

figure
plot(x(1:N),ruy_norm)
hold all
grid on
plot(x(1:N),y1(1:N),x(1:N),u(1:N))
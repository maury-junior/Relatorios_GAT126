load('ensaio_prbs.txt')

t = ensaio_prbs(:,1);
u = ensaio_prbs(:,2);
y = ensaio_prbs(:,3);

% figure('color',[1 1 1])
% plot(t,u)


N = length(u)/2;

ry_norm = zeros(N,1);
ry2_norm = zeros(N,1);
for k=1:N
    ry_T_norm = (y(1:N)-mean(y))'*(y(1+k:N+k)-mean(y));
    ry_norm(k) = ry_T_norm/(2*N+1);
    
    ry2_T_norm = (y(1:N).^2-mean(y.^2))'*(y(1+k:N+k).^2-mean(y.^2));
    ry2_norm(k) = ry2_T_norm/(2*N+1);
    
    
end

figure
plot(t(1:N),ry_norm)

ruu = zeros(N,1);
for k=1:N
    ruu_T = u(1:N)'*u(1+k:N+k);
    ruu(k) = ruu_T/(2*N+1);
end

ruy = zeros(N,1);
for k=1:N
    ruy_T = u(1:N)'*y(1+k:N+k);
    ruy(k) = ruy_T/(2*N+1);
end

h_FAC = ruy./ruu;
h_FAC_f = fft(ruy)./fft(ruu);

freq = 2*pi*1/N*(0:N/2);

figure
semilogx(freq,20*log10(abs(h_FAC_f(1:length(freq)))))

% function ru = AC(u,k)
% 
% if size(u,1) == 1 && size(u,1) > 1
%     u = u';
% elseif size(u,1) == 1 && size(u,1) == 1
%     error('u is a scalar.')
% elseif size(u,1) > 1 && size(u,2) > 1
%     error('u must be a vector, not a matrix.')
% end
% 
% ruu = 
% 
% end
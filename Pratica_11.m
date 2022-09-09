load('ensaio_prbs.txt')

t_dados = ensaio_prbs(:,1);
u_dados = ensaio_prbs(:,2);
y_dados = ensaio_prbs(:,3);

% figure('color',[1 1 1])
% plot(t,u)


N = length(u_dados)/2;

ry_norm = zeros(N,1);
ry2_norm = zeros(N,1);
for k=1:N
    ry_T_norm = (y_dados(1:N)-mean(y_dados))'*(y_dados(1+k:N+k)-mean(y_dados));
    ry_norm(k) = ry_T_norm/(2*N+1);
    
    ry2_T_norm = (y_dados(1:N).^2-mean(y_dados.^2))'*(y_dados(1+k:N+k).^2-mean(y_dados.^2));
    ry2_norm(k) = ry2_T_norm/(2*N+1);
end

figure('color',[1 1 1])
subplot(211)
plot(t_dados(1:N),ry_norm)
subplot(212)
plot(t_dados(1:N),ry2_norm)

%%
d = 40;

t_d = t_dados(1:d:length(t_dados));
u_d = u_dados(1:d:length(t_dados));
y_d = y_dados(1:d:length(t_dados));

N_d = length(t_d)/2;
ry_d_norm = zeros(N_d,1);
ry2_d_norm = zeros(N_d,1);
for k=1:N_d
    ry_T_d_norm = (y_d(1:N_d)-mean(y_d))'*(y_d(1+k:N_d+k)-mean(y_d));
    ry_d_norm(k) = ry_T_d_norm/(2*N_d+1);
    
    ry2_T_d_norm = (y_d(1:N_d).^2-mean(y_d.^2))'*(y_d(1+k:N_d+k).^2-mean(y_d.^2));
    ry2_d_norm(k) = ry2_T_d_norm/(2*N_d+1);
end

figure('color',[1 1 1])
subplot(211)
plot(t_dados(1:N_d),ry_d_norm)
subplot(212)
plot(t_dados(1:N_d),ry2_d_norm)

%%
t = t_d;
u = u_d;
y = y_d;

N = floor(length(t)/2);

ruu_norm = zeros(N,1);
for k=1:N
    ruu_norm_T = (u(1:N)'-mean(u))*(u(1+k:N+k)-mean(u));
    ruu_norm(k) = ruu_norm_T/(2*N+1);
end

ruy_norm = zeros(N,1);
for k=1:N
    ruy_norm_T = (u(1:N)'-mean(u))*(y(1+k:N+k)-mean(y));
    ruy_norm(k) = ruy_norm_T/(2*N+1);
end

figure
plot(t(1:N),ruy_norm)

h_FAC = ruy_norm./ruu_norm;
h_FAC_f = fft(ruy_norm)./fft(ruu_norm);

freq = 2*pi*1/N*(0:N/2);

figure
semilogx(freq,20*log10(abs(h_FAC_f(1:length(freq)))))

%%


ruu_norm = zeros(N,1);
for k=1:N
    ruu_norm_T = (u_dados(1:N)'-mean(u_dados))*(u_dados(1+k:N+k)-mean(u_dados));
    ruu_norm(k) = ruu_norm_T/(2*N+1);
end

ruy_norm = zeros(N,1);
for k=1:N
    ruy_norm_T = (u_dados(1:N)'-mean(u_dados))*(y_dados(1+k:N+k)-mean(y_dados));
    ruy_norm(k) = ruy_norm_T/(2*N+1);
end


figure
plot(t_dados(1:1500),ruy_norm)

figure
plot(t_dados,u_dados,'b',t_dados,y_dados,'r')

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
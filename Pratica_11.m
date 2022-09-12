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
plot(1:N,ry_norm)
subplot(212)
plot(1:N,ry2_norm)

ruy = zeros(N,1);
for k=1:N
    ruy_T = (u_dados(1:N)'-mean(u_dados))*(y_dados(1+k:N+k)-mean(y_dados));
    ruy(k) = ruy_T/(2*N+1);
end

figure('color',[1 1 1])
plot(ruy)
grid on
% atraso = 4
%%
d = 6; % Dizimização

t_d = t_dados(1:d:length(t_dados));
u_d = u_dados(1:d:length(t_dados));
y_d = y_dados(1:d:length(t_dados));

N_d = floor(length(t_d)/2);
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
plot(1:N_d,ry_d_norm)
subplot(212)
plot(1:N_d,ry2_d_norm)

%% Método MQ
t = t_d;
u = u_d;
y = y_d;

tau_d = 5; % atraso

% y = a y(k-1) + b x(k-4) + c x(k-5)
N_d = length(t);
phi1 = [y(tau_d+2:N_d-1) u(2:N_d-tau_d-1) u(1:N_d-tau_d-2)];
theta1 = pinv(phi1)*y(tau_d+3:N_d);

phi2 = [y(tau_d+3:N_d-1) y(tau_d+2:N_d-2) u(3:N_d-tau_d-1) u(2:N_d-tau_d-2) u(1:N_d-tau_d-3)];
theta2 = pinv(phi2)*y(tau_d+4:N_d);

% Simulação livre
% 1 atraso
Y1 = zeros(length(y),1);
Y1(1:tau_d+3) = y(1:tau_d+3);
for k=tau_d+4:N_d
    Y1(k) = [Y1(k-1) u(k-4) u(k-5)]*theta1;
end

% 2 atrasos
Y2 = zeros(length(y),1);
Y2(1:tau_d+4) = y(1:tau_d+4);
for k=tau_d+5:N_d
    Y2(k) = [Y2(k-2) Y2(k-1) u(k-1) u(k-2) u(k-3)]*theta2;
end

figure('color',[1 1 1])
plot(y)
hold on
plot(Y1,'k')
plot(Y2,'r')
legend('data','Um atraso','Dois atrasos')

MSE = [sum((y-Y1).^2) sum((y-Y2).^2)]

%% Método FACFCC
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

% figure
% plot(1:N,ruy_norm)
Ruu = zeros(125,125);

for i=1:125
   Ruu(i,:) = ruu_norm(124+i:-1:i,1)';
end

h_FAC = inv(Ruu)*ruy_norm(125:249);

figure('color',[1 1 1])
plot(h_FAC)

u_deg = ones(125,1);
U_step = u_deg; % criacao da matrix de entrada U
for i=1:length(u_deg)-1
    U_step = [U_step [zeros(i,1); u_deg(1:length(u_deg)-i)]];
end;

Y_deg = U_step*h_FAC;

h_FAC_f = fft(ruy_norm)./fft(ruu_norm);

freq = 2*pi*1/N*(0:N/2);

figure
semilogx(freq,20*log10(abs(h_FAC_f(1:length(freq)))))
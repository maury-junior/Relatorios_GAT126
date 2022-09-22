%% Parte 1
dados = load('saida_degrau.txt');
t = dados(:,1);
v = dados(:,2);
T = dados(:,3);

y0 = mean(T(1:20));
yinf = mean(T(end-20:end));
du = mean(v);
K = (yinf-y0)/du;

s = tf('s');
num = [K];

% Modelo de primeira ordem
%a - Metodo das integrais
% usando trapz
v_norm = (v-0)./(max(v) - 0);
T_norm = (T-y0)./(yinf - y0);
theta_tau = trapz(t,v_norm-T_norm);
tau_1_Nish = exp(1)*trapz(t(t<theta_tau),T_norm(t<theta_tau));
t_d_Nish = theta_tau-tau_1_Nish;
den_Nish = [tau_1_Nish 1];
H_Nish = tf(num,den_Nish);
H_Nish = H_Nish*exp(-t_d_Nish*s)

%b - Método de Mollenkamp
y1 = y0 + 0.15*(yinf-y0);
y2 = y0 + 0.45*(yinf-y0);
y3 = y0 + 0.75*(yinf-y0);

for i=1:length(T)-10
    T_f(i) = mean(T(i:i+10));
end

ind = find(T>y1);
t1 = t(ind(1)); % 0.15 y(inf)
ind = find(T>y2);
t2 = t(ind(1)); % 0.45 y(inf)
ind = find(T>y3);
t3 = t(ind(1));% 0.75 y(inf)

H_Mol = SOPDT_Mollenkamp(K,t1,t2,t3)

% figure('color',[1 1 1])
% plot(t,T,'r')
% hold all
% step(H_Nish,t)
% step(H_Mol,t)

dados = load('saida_prbs_parte1.txt');
t = dados(:,1);
v = dados(:,2);
T = dados(:,3);

T_Nish = lsim(H_Nish,v,t);
T_Mol = lsim(H_Mol,v,t);
% figure('color',[1 1 1])
% plot(t,T);
% hold all
% plot(t,T_Nish)
% plot(t,T_Mol)
% legend('dados','Nish','Mollenkamp')
MSE_Nish = sum((T-T_Nish).^2)
MSE_Mol = sum((T-T_Mol).^2)

%% Parte 2

dados = load('saida_prbs_parte2.txt');
t = dados(:,1);
v = dados(:,2);
T = dados(:,3);

% autocorr(v)
[ruy, lags] = crosscorr(v,T,1499);
n_d = lags(find(ruy==max(ruy)))
tau_d = n_d*(t(2)-t(1));

figure
crosscorr(v,T,1499)

[ry,lagsy] = autocorr(T,1499);
indy = find(islocalmin(ry));
lag = lagsy(indy(1));
% [ry2,lagsy2] = autocorr(T.^2,1499);
% indy2 = find(islocalmin(ry2));
% lag = min(lagsy(indy(1), lagsy2(indy2(1)))

d = floor(lag/10)


t_d = decimate(t,d); 
v_d = decimate(v,d);
T_d = decimate(T,d);
% t_d = t(1:d:end);
% v_d = v(1:d:end);
% T_d = T(1:d:end);
T_amostragem = t_d(2)-t_d(1)
% Identificacao do sistema pelo metodo da convolucao
U=v_d; % criacao da matrix de entrada U
for i=1:length(v_d)-1
%     U = [U [zeros(i,1); v_d(1:length(v_d)-i)]];
    U = [U [v_d(length(v_d)-i:end-1); v_d(1:length(v_d)-i)]];
end

H=inv(U)*T_d;

figure
plot(t_d,H)
grid on
hold all

% FACFCC
figure
autocorr(T_d,length(T_d)-1)
[ruu, lagsu] = autocorr(v_d,length(v_d)-1);
[ruy, lagsuy] = crosscorr(v_d,T_d,length(v_d)-1);

H = ruy(length(v_d):end)/var(ruu);

H_f = fft(ruy(length(v_d):end))./fft(ruu);
freq = 2*pi/length(ruu)*(1:length(ruu)/2);

figure
subplot(211)
semilogx(2*pi*freq,20*log10(abs(H_f(1:length(freq)))),'k--');
title('(a)');
ylabel('|H(j\omega)| (dB)');
subplot(212)
semilogx(2*pi*freq,angle(H_f(1:length(freq)))*180/pi,'k--');
title('(b)');
ylabel('fase de H(j\omega) (graus)');
xlabel('\omega (rad/s)')

%% 
dados = load('saida_prbs_parte2.txt');
t = dados(:,1);
v = dados(:,2);
T = dados(:,3);

% autocorr(v)
[ruy, lags] = crosscorr(v,T,1499);
n_d = lags(find(ruy==max(ruy)))
tau_d = n_d*(t(2)-t(1));

[ry,lagsy] = autocorr(T,1499);
indy = find(islocalmin(ry));
lag = lagsy(indy(1));
% [ry2,lagsy2] = autocorr(T.^2,1499);
% indy2 = find(islocalmin(ry2));
% lag = min(lagsy(indy(1), lagsy2(indy2(1)))

d = floor(lag/10)

t_d = t(1:d:end);
v_d = v(1:d:end);
T_d = T(1:d:end);

% Modelo 1 atraso -> Ya(k) = [T_d(k-1) v_d(k-1-n_d)]* theta1
phi1 = [T_d(n_d:end-1) v_d(1:end-n_d)];
Y1 = T_d(n_d+1:end);
theta1 = pinv(phi1)*Y1;

T1(1:n_d) = T_d(1:n_d);
for k=n_d+1:length(T_d)
    T1(k) = theta1(1)*T_d(k-1) + theta1(2)*v_d(k-n_d);
end

% Modelo 2 atrasos -> Ya(k) = [T_d(k-1) T_d(k-2) v_d(k-1-n_d) v_d(k-2-n_d)]* theta1
phi2 = [T_d(n_d+1:end-1) T_d(n_d:end-2) v_d(2:end-n_d) v_d(1:end-n_d-1)];
Y2 = T_d(n_d+2:end);
theta2 = pinv(phi2)*Y2;

T2(1:n_d+1) = T_d(1:n_d+1);
for k=n_d+2:length(T_d)
    T2(k) = theta2(1)*T_d(k-1) + theta2(2)*T_d(k-2) + theta2(3)*v_d(k-n_d) + theta2(4)*v_d(k-n_d-1);
end

figure('color',[1 1 1])
plot(t_d,T_d,'b')
hold all
plot(t_d,T1,'k')
plot(t_d,T2,'r')

MSE_1 = mean((T1'-T_d).^2)
MSE_2 = mean((T2'-T_d).^2)

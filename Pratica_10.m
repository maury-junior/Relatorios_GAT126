%%  Parte 1 - Sistema sobreamortecido
s = tf('s');

dadosordem2sob = load('dadosordem2sobre.txt');
t_sob = dadosordem2sob(:,1);
v_sob = dadosordem2sob(:,2);
T_sob = dadosordem2sob(:,3);

Dx_sob = max(v_sob)-0;
Dy_sob = max(T_sob)-min(T_sob);
  
K_sob = Dy_sob/Dx_sob;

% Mollenkamp

y_inf_sob = 10;

t1 = 42;    % 0.15 y(inf)
t2 = 83.5;  % 0.45 y(inf)
t3 = 144.5; % 0.75 y(inf)

H_Mol_sob = SOPDT_Mollenkamp(K_sob,t1,t2,t3);

% Caso sobreamortecido
u_norm_sob = v_sob./v_sob;
y_norm_sob = (T_sob-min(T_sob))./(max(T_sob)-min(T_sob));

figure('color',[1 1 1])
plot(t_sob,y_norm_sob);
grid on

Mi = 0.4/50;
tm = 150;
sys_type = 1; % sistema sobreamortecido
H_Sun_sob = SOPDT_Sundaresan(K_sob,t_sob,u_norm_sob,y_norm_sob,sys_type,Mi,tm);
% eta = 0.33

figure('color',[1 1 1])
plot(t_sob,T_sob,'r')
grid on
hold on
step(H_Mol_sob,t_sob)
step(H_Sun_sob,t_sob)
%% Parte 2 - Sistema subamortecido
dadosordem2sub = load('dadosordem2sub.txt');
t_sub = dadosordem2sub(:,1);
v_sub = dadosordem2sub(:,2);
T_sub = dadosordem2sub(:,3);

Dx_sub = max(v_sub)-0;
Dy_sub = 5-min(T_sub);
K_sub = Dy_sub/Dx_sub;

td_sub = 3.02;
tp_sub = 4.68-td_sub;
M_p = (6.86-5)/5;

figure('color',[1 1 1])
plot(t_sub,T_sub,'r')
hold on 
grid on

H_Mp = SOPDT_Mp_equations(K_sub,tp_sub,M_p,td_sub);
step(H_Mp)
%% Mollenkamp

% y_inf_sob = 10;
% 
% t1 = 42;    % 0.15 y(inf)
% t2 = 83.5;  % 0.45 y(inf)
% t3 = 144.5; % 0.75 y(inf)
% 
% H_Mol_sob = SOPDT_Mollenkamp(K_sob,t1,t2,t3);

y_inf_sub = 5;
t1 = 3.29; % 0.15 y(inf)
t2 = 3.56; % 0.45 y(inf)
t3 = 3.78; % 0.75 y(inf)

H_Mol_sub = SOPDT_Mollenkamp(K_sub,t1,t2,t3);

%% Philipp e Parr
% Number of visible cycles
N = 2;
TN = 11.29-4.68;
H_PP = SOPDT_Phillip_Parr(K_sub,N,TN,td_sub);
step(H_PP)
%% Sundaresan
% % Caso sobreamortecido
% u_norm_sob = v_sob./v_sob;
% y_norm_sob = (T_sob-min(T_sob))./(max(T_sob)-min(T_sob));
% 
% figure('color',[1 1 1])
% plot(t_sob,y_norm_sob);
% grid on
% 
% Mi = 0.4/50;
% tm = 150;
% sys_type = 1; % sistema sobreamortecido
% H_Sun_sob = SOPDT_Sundaresan(K_sob,t_sob,u_norm_sob,y_norm_sob,sys_type,Mi,tm);

% Caso subamortecido
u_norm_sub = v_sub./v_sub;
y_norm_sub = (T_sub-min(T_sub))./(5-min(T_sub));

figure('color',[1 1 1])
plot(t_sub,y_norm_sub);
grid on

Mi = (0.996-0.204)/(3.98-3.35);
tm = 3.98;
sys_type = 2; % sistema subamortecido
H_Sun_sub = SOPDT_Sundaresan(K_sub,t_sub,u_norm_sub,y_norm_sub,sys_type,Mi,tm);

%% Parte 3
dadosordem2rand = load('dadosordem2aleatorio.txt');
t_rand = dadosordem2rand(:,1);
v_rand = dadosordem2rand(:,2);
T_rand = dadosordem2rand(:,3);


% Identificacao pelo metodo da convolucao
plot(t_rand,T_rand)
grid on
hold on

% Identificacao do sistema pelo metodo da convolucao
U=v_rand; % criacao da matrix de entrada U
for i=1:length(v_rand)-1
    U = [U [zeros(i,1); v_rand(1:length(v_rand)-i)]];
end;

% solucao do sistema Y=UX para U (U = Y/X)
% H1=U\T_rand;
H=inv(U)*T_rand;

% figure
% plot(H1)
figure
plot(t_rand,H)
axis([0 5 -2.6 2.6])
grid on
hold all

u = ones(2000,1);
U_step = u; % criacao da matrix de entrada U
for i=1:length(u)-1
    U_step = [U_step [zeros(i,1); u(1:length(u)-i)]];
end;

Y_st = U_step*H;

plot(Y_st)

% Identificacao da resposta em frequencia
v_rand_d = v_rand(1:length(v_rand)-1) - v_rand(2:length(v_rand));
T_rand_d = T_rand(1:length(T_rand)-1) - T_rand(2:length(T_rand));

v_d_fft = fft(v_rand_d);
T_d_fft = fft(T_rand_d);
H = T_d_fft./v_d_fft;

% vetor de frequencias
freq=1/(length(T_rand))*(0:length(T_rand)/2);

figure
subplot(211)
semilogx(2*pi*freq,20*log10(abs(H(1:length(freq)))),'k--');
title('(a)');
ylabel('|H(j\omega)| (dB)');
subplot(212)
semilogx(2*pi*freq,angle(H(1:length(freq)))*180/pi,'k--');
title('(b)');
ylabel('fase de H(j\omega) (graus)');
xlabel('\omega (rad/s)')
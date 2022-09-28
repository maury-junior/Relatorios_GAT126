%% Parte 3 - Resposta ao impulso e resposta em frequencia
dadosordem2rand = load('dadosordem2aleatorio.txt');
t_rand = dadosordem2rand(:,1);
v_rand = dadosordem2rand(:,2);
T_rand = dadosordem2rand(:,3);

% Identificacao pelo metodo da convolucao
plot(t_rand,T_rand)
grid on
hold on

for i=1:length(v_rand)-8
    v_f(i) = mean(v_rand(i:i+8));
    T_f(i) = mean(T_rand(i:i+8));
end
t_f = t_rand(1:length(v_rand)-8);
v_f = v_f';
T_f = T_f';

plot(t_f,T_f,'r')

% figure
% plot(t_f,v_f)

[ruy lagsuy]=crosscorr(v_f,T_f,1991);
[ry lagsy]=autocorr(T_f,1998);
% [ruy lagsuy]=xcorr(v_rand,T_rand);
% [ry lagsy]=xcorr(T_rand);

% figure('color',[1 1 1])
% plot(lagsuy,ruy)
% figure('color',[1 1 1])
% plot(lagsy,ry)

U=v_f; % criacao da matrix de entrada U
for i=1:length(v_f)-1
    U = [U [zeros(i,1); v_f(1:length(v_f)-i)]];
end;

% solucao do sistema Y=UX para U (U = Y/X)
% H1=U\T_rand;
H=inv(U)*T_f;

% figure
% plot(H1)
figure
plot(t_f,H)
axis([0 5 -2.6 2.6])
grid on
hold all

H_r = ruy(1980:end-12)/var(v_f);
figure
plot(t_f,H_r)


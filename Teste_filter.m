%% Parte 3 - Resposta ao impulso e resposta em frequencia
dadosordem2rand = load('dadosordem2aleatorio.txt');
t_rand = dadosordem2rand(:,1);
v_rand = dadosordem2rand(:,2);
T_rand = dadosordem2rand(:,3);

% Identificacao pelo metodo da convolucao
plot(t_rand,T_rand)
grid on
hold on

[ruy lagsuy]=crosscorr(v_rand,T_rand);
[ry lagsy]=autocorr(T_rand);
% [ruy lagsuy]=xcorr(v_rand,T_rand);
% [ry lagsy]=xcorr(T_rand);

figure('color',[1 1 1])
plot(lagsuy,ruy)
figure('color',[1 1 1])
plot(lagsy,ry)


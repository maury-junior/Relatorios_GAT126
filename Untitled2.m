dadosordem2rand = load('dadosordem2aleatorio.txt');
t_rand = dadosordem2rand(:,1);
v_rand = dadosordem2rand(:,2);
T_rand = dadosordem2rand(:,3);

plot(t_rand,T_rand)
grid on
hold on

% Filtro
for i = 1:length(T_rand)-8
    T_rand_f(i) = mean(T_rand(i:i+8));
    v_rand_f(i) = mean(v_rand(i:i+8));
end
t_rand_f = t_rand(1:length(T_rand)-8)';
% v_rand_f = v_rand(1:length(T_rand)-8)';

plot(t_rand_f,T_rand_f,'r')

% Identificacao do sistema pelo metodo da convolucao
U=v_rand_f'; % criacao da matrix de entrada U
for i=1:length(v_rand_f)-1
    U = [U [zeros(i,1); v_rand_f(1:length(v_rand_f)-i)']];
end

% solucao do sistema Y=UX para U (U = Y/X)
H1=U\T_rand_f';
H2=inv(U)*T_rand_f';

figure
plot(t_rand_f,H1)
grid on
axis([0 10 -5 5])
figure
plot(t_rand_f,H2)

U2=v_rand; % criacao da matrix de entrada U
for i=1:length(v_rand)-1
    U2 = [U2 [zeros(i,1); v_rand(1:length(v_rand)-i)]];
end

H = U2\T_rand;

figure
plot(t_rand,H)
grid on
axis([0 10 -5 5])
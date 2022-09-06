dadosordem2sub = load('dadosordem2sub.txt');
t_rand = dadosordem2sub(:,1);
v_rand = dadosordem2sub(:,2);
T_rand = dadosordem2sub(:,3);

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
hold all

u = ones(1000,1);
U_step = u; % criacao da matrix de entrada U
for i=1:length(u)-1
    U_step = [U_step [zeros(i,1); u(1:length(u)-i)]];
end;

Y_st = U_step*H;

plot(t_rand,Y_st)
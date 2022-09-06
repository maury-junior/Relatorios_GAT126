dados = load('dadosordem2aleatorio.txt');

%% Convolução

for i=1:1000
    A(i,:)=dados(999+i:-1:i,2);
end

h = inv(A)*dados(1000:1999,3);

yi = dimpulse([1 0.5],[1 -1.5 0.7],length(h));

figure
plot(h,'o-')
hold on
plot(yi,'r')
 
%% Parte 3

dadosordem2rand = load('dadosordem2aleatorio.txt');
t_rand = dadosordem2rand(:,1);
v_rand = dadosordem2rand(:,2);
T_rand = dadosordem2rand(:,3);


% Identificacao pelo metodo da convolucao
figure('color',[1 1 1])
plot(t_rand(1000:1999),T_rand(1000:1999),'r')
grid on
hold on

% Identificacao do sistema pelo metodo da convolucao
% U=v_rand; % criacao da matrix de entrada U
for i=1:1000
    U(i,:) = v_rand(999+i:-1:i);
end;

% solucao do sistema Y=UX para U (U = Y/X)
% H1=U\T_rand;
H=inv(U)*T_rand(1000:1999);

% figure
% plot(H1)
figure
plot(t_rand(1000:1999),H)
hold all

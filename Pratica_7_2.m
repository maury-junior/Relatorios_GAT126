%% b)
m1=1;
clear
clc

m1=1;
m2=2;
k1=100;
k2=20;
k3=20;
b=12;

% X1/U
num1 = [m2 b k2+k3];
den = [m1*m2, b*(m1+m2), m1*(k2+k3)+m2*(k1+k2), b*(k1+k3) (k1*k2+k1*k3+k2*k3)];
G1 = tf(num1,den)

% X2/U
num2 = [b k2];
G2 = tf(num2,den)

%% Degrau
figure
step(G1)
grid on
hold on
step(G2)
legend('G1','G2')

%%
der = tf('s')
V1 = G1*der
V2 = G2*der

A1 =G1*der^2
A2 = G2*der^2
    
step(G1); 
step(V1); 
step(A1); 

%legend('G1','V1','A1')
% ylabel({'Posicao [m]','Velocidade [m/2]','Aceleracao [m/s^2]'})

%% Impulso
impulse(G1)
grid on
hold on
impulse(G2)
legend('G1','G2')

%% Senoidal
t = 0:0.01:5;
u = sin(pi*t);
lsim(G1,u,t)
grid on
hold on
lsim(G2,u,t)
legend('G1','G2')

%% Espaço de estados
[A1,B1,C1,D1] = tf2ss(num1,den)
G1ss = ss(A1,B1,C1,D1)

[A2,B2,C2,D2]= tf2ss(num2,den)
G2ss = ss(A2,B2,C2,D2)

figure
t = 0:0.01:5;
u = sin(pi*t);
lsim(G1ss,u,t)
grid on
hold on
lsim(G2ss,u,t)
legend('G1ss','G2ss')



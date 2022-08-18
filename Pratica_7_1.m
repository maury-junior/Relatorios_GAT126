%% b)
k = 1;
b = 1;
m = 1;

num1 = [b k];
den1 = [m b k];

G1 = tf(num1,den1)
step(G1)
hold on
grid on

%% c)
k2 = 4*k;
num2 = [b k2];
den2 = [m b k2];
G2 = tf(num2,den2)
step(G2)
%% d)
b3 = 4*b;
num3 = [b3 k];
den3 = [m b3 k];
G3 = tf(num3,den3)
step(G3)

%% e)
figure
bode(G1)
grid on
hold on
bode(G2)
bode(G3)
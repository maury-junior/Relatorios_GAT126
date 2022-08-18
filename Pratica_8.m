%% 1
% b)
C = 2;
R = 10;

H_Qi = tf([R],[R*C 1])
Qo_Qi = tf([1],[R*C 1])

step(H_Qi)
hold on
step(Qo_Qi)

legend('h(t)','Qo(t)')
% e)
figure
bode(H_Qi)
hold on
bode(Qo_Qi)
grid on

%% 2
C_1 = 1.5;
C_2 = 2;
R_1 = 1;
R_2 = 5;
G11 = tf([(R_1^2*R_2*C_2) 0 (R_1^2+R_1*R_2)],[(R_1^2*C_1*R_2*C_2) (R_1*R_2*C_2+R_1^2*C_1+R_1*R_2*C_1) R_1])

G12 = tf([R_1*R_2],[(R_1^2*C_1*R_2*C_2) (R_1*R_2*C_2+R_1^2*C_1+R_1*R_2*C_1) R_1])

G21 = tf([R_2],[(R_1*C_1*R_2*C_2) (R_2*C_2+R_1*C_1+R_2*C_1) 1])

% G22 = tf([R_1*R_2],[(R_1^2*C_1*R_2*C_2) (R_1*R_2*C_2+R_1^2*C_1+R_1*R_2*C_1) R_1])
G22 = tf([R_1*R_2*C_1 R_2],[(R_1*C_1*R_2*C_2) (R_2*C_2+R_1*C_1+R_2*C_1) 1])

figure
t = 0:0.5:200;
[h11,t]=step(G11,t);
opt = stepDataOptions('StepAmplitude',3);
[h12,t]=step(G12,t,opt);
h1 = h11+h12;
plot(t,h1)
grid on
hold on

[h21,t]=step(G21,t);
[h22,t]=step(G22,t,opt);
h2 = h21+h22;
plot(t,h2)

legend('H1','H2')
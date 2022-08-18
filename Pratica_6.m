%% Pratica 6
%% 1)
R = 100;
C = 1e-3;
G1 = tf([1],[R*C 1]) % Vc/Ve
step(G1)
ylabel('V_c(t)')
G1_2R = tf([1],[2*R*C 1]) % Vc/Ve para dobro de R

G2 = tf([C 0],[R*C 1]) % i/Ve
figure
step(G2)
ylabel('i(t)')

%% 2)
R = 100;
L = 1e-3;

G3 = tf([1],[L R]) % I/Ve
figure
step(G3)
ylabel('i(t)')

G4 = tf([L 0],[L R]) % V_L/Ve
figure
step(G4)
ylabel('V_L(t)')

%% 3)
R = 100;
L = 1e-3;
C = 1e-3;

G5 = tf([1],[L*C R*C 1]) % V_C/Ve
figure
step(G5)
ylabel('V_C(t)')

G6 = tf([C 0],[L*C R*C 1]) % I/Ve
figure
step(G6)
ylabel('i(t)')

%% 4)
figure
bode(G1)
hold on
bode(G4)
bode(G5)
legend('V_c RC','V_L RL','V_c RLC')

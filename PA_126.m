% Parte 1
H1 = tf([1],[5 1])*exp(-5*s);
[y_P1,t] = step(H1);
v = 1e-2*(randn(length(y_P1),1)-0.5);
figure('color',[1 1 1])
plot(t,y_P1+v)
grid on

% Parte 2
zeta = 0.15; wn = 3;
H2 = tf([wn^2],[1 2*zeta*wn wn^2])*exp(-2*s);

[y_P2,t] = step(H2);
v = 1e-2*(randn(length(y_P2),1)-0.5);
figure('color',[1 1 1])
plot(t,y_P2+v)
grid on
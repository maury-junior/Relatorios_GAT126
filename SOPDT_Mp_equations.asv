function H = SOPDT_Mp_equations(K,tp,M_p,td)

s = tf('s');
zeta = -log(M_p)/sqrt(log(M_p)^2+pi^2);

wn = pi/(tp*sqrt(1-zeta^2));

H = tf(K*wn^2,[1 2*zeta*wn wn^2])*exp(-td*s);
function H = SOPDT_Phillip_Parr(K,N,TN,td)
% K - ganho
% N - Numero de ciclos visiveis
% TN - Tempo para os N ciclos
s = tf('s');
zeta = 0.6/N;
T = TN/N; % periodo medio
wn = 2*pi/T;

teste_Ph_Parr = sqrt(1-zeta^2);
if teste_Ph_Parr < 0.97
    disp('Philip and Parr method not appropriate' )
end

H = tf(K*wn^2,[1 2*zeta*wn wn^2])*exp(-td*s);


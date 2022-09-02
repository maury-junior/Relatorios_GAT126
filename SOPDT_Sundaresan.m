function H = SOPDT_Sundaresan(K,t,u_norm,y_norm,sys_type,Mi,tm)

% u_norm   - entrada normalizada
% y_norm   - saida normalizada

% sys_type - type of system, 1- sobreamortecido, 2- subamortecido

% Mi       - Inclinacao da tangente no pto de inflexao de y(t)
% tm       - Intersecao da tangente com o valor em reg permanente
s = tf('s')

if sys_type == 1 % caso sobreamortecido
    
    m1 = trapz(t,u_norm-y_norm);
    lambda = (tm-m1)*Mi;
    sprintf('lambda = %d',lambda)
    prompt = 'Enter eta: ';
    eta = input(prompt); % Obtido a partir de grafico (fig 3.4)
    tau_1 = eta^(eta/(1-eta))/Mi;
    tau_2 = eta^(1/(1-eta))/Mi;
    tau_d = m1 - tau_1 - tau_2;
    
    wn = sqrt(1/(tau_1*tau_2));
    zeta = (tau_1+tau_2)/(2*sqrt(tau_1*tau_2));
    
elseif sys_type == 2 % caso subamortecido
    m1 = trapz(t,u_norm-y_norm);
    lambda = (tm-m1)*Mi;
    sprintf('lambda = %d',lambda)
    prompt = 'Enter zeta: ';
    zeta = input(prompt); % Obtido a partir de grafico (fig 3.6)
    wn = acos(zeta)/sqrt(1-zeta^2) * 1/(tm-m1);
    tau_d = m1 - 2*zeta/wn;
end

H = tf(K*wn^2,[1 2*zeta*wn wn^2])*exp(-tau_d*s);

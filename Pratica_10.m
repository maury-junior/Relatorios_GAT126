Dx = 1;
Dy = 1;

M_p = 0.3;

K = Dy/Dx;


%% Mollenkamp

t1 % 0.15 y(inf)
t2 % 0.45 y(inf)
t3 % 0.75 y(inf)

x = (t2-t1)/(t3-t1);
zeta = (0.0805-5.547*(0.475-x)^2)/(x-0.356);
if zeta <1
   f1=0.708*2.811^zeta ;
elseif zeta >=1
    f1=2.6*zeta-0.6;
end

omega_n = f1/(t3-t1);

f2=0.922*1.66^zeta;
theta = t2-f2/omega_n;

if zeta>=1
    tau_1 = (zeta+sqrt(zeta^2-1))/omega_n;
    tau_2 = (zeta-sqrt(zeta^2-1))/omega_n;
end

%% Philipp e Parr
% Number of visible cycles
N = 4;
omega_n = 1;

zeta = 0.6/N;

teste_Ph_Parr = sqrt(1-zeta^2);
if teste_Ph_Parr > 0.2
    disp('Philip and Parr method not appropriate' )
end

%% Sundaresan
% K = Dy/Dx;
u_norm % entrada normalizada
y_norm % saida normalizada

if zeta >= 1
    Mi % Inclinacao da tangente no pto de inflexao de y(t)
    tm % Intersecao da tangente com o valor em reg permanente
    m1 = trapz(t,u_norm-y_norm);
    lambda = (tm-m1)*Mi;
    eta % Obtido a partir de grafico (fig 3.4)
    tau_1 = eta^(eta/(1-eta))/Mi;
    tau_2 = eta^(1/(1-eta))/Mi;
    tau_d = m1 - tau_1 - tau_2;
    
    tau_1_Nish = exp(1)*trapz(t(t<theta_tau),T_norm(t<theta_tau));

elseif zeta < 1
    tm
    Mi
    m1 = trapz(t,u_norm-y_norm);
    lambda = (tm-m1)*Mi;
    zeta % Graficamente (fig 3.6)
    omega_n = acos(zeta)/sqrt(1-zeta^2) * 1/(tm-m1);
    td = m1 - 2*zeta/omega_n;
end
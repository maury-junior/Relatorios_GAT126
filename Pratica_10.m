
dados = load('dadosordem2sobre.txt');
t = dados(:,1);
v = dados(:,2);
T = dados(:,3);

figure('color',[1 1 1])
plot(t,T)
grid on
hold on


Dx = max(v) - 0;
Dy = max(T) - min(T);

K = Dy/Dx;

M_p = 0.3;

%% Mollenkamp

y_15 = min(T) + 0.15*Dy;
y_45 = min(T) + 0.45*Dy;
y_75 = min(T) + 0.75*Dy;

plot(t,y_15*ones(1,length(t)))
plot(t,y_45*ones(1,length(t)))
plot(t,y_75*ones(1,length(t)))

%

t1 =  42;   % 0.15 y(inf)
t2 =  83;   % 0.45 y(inf)
t3 = 143.5; % 0.75 y(inf)

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

H_Mol = tf([K*omega_n^2],[1 2*zeta*omega_n omega_n^2]);
H_Mol_d = H_Mol*exp(-11*s);

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
% Pratica_9_Metodos_deterministicos_1_ordem
dados = load('dadosordem1.txt')
t = dados(:,1);
v = dados(:,2);
T = dados(:,3);

figure('color',[1 1 1])
plotyy(t,v,t,T)
grid on
legend('V(t)','T(t)')

K_approx = (8-2)/2;
t_d_approx = 5;
y_632_approx = 2+0.632*(8-2);
t_632_approx = 15;
s = tf('s');
num = [K_approx];

%a - Metodo das integrais
%usando trapz
v_norm = (v-0)./(max(v) - 0)
T_norm = (T-min(T))./(max(T) - min(T))
theta_tau = trapz(t,v_norm-T_norm)
tau_1_Nish = exp(1)*trapz(t(t<theta_tau),T_norm(t<theta_tau));
t_d_Nish = theta_tau-tau_1_Nish;
den_Nish = [tau_1_Nish 1];
H_Nish = tf(num,den_Nish);
H_Nish = H_Nish*exp(-t_d_Nish*s);
% step(H_Nish)
%H_Nish = tf(num,den_Nish,'InputDelay',t_d_Nish)

%b - Metodo de uma constante de tempo (63,2%)
tau_1_1const_temp = 15-t_d_approx;
t_d_1const_temp = t_d_approx;
den_1const = [tau_1_1const_temp 1];
H_1const = tf(num,den_1const)*exp(-t_d_1const_temp*s);

%c - Metodo de quatro constantes de tempo (98%)
y_4t_approx = 2+0.98*(8-2);
t_4t_approx = 42;
tau_4_approx = t_4t_approx - t_d_approx;
tau_1_4const_temp = tau_4_approx/4;
t_d_4const_temp = t_d_approx;
den_4const = [tau_1_4const_temp 1];
H_4const = tf(num,den_4const)*exp(-t_d_4const_temp*s);

%d - Metodo de Ziegler-Nichols;
t_3_ZN = 16.7;
tau_1_ZN = t_3_ZN-t_d_approx;
t_d_ZN = t_d_approx;
den_ZN = [tau_1_ZN 1];
H_ZN = tf(num, den_ZN)*exp(-t_d_ZN*s);

%e - Metodo de Smith;
y_283_approx = 2 + 0.283*(8-2);
t_1_Sm = 8.3;
t_2_Sm = 15;
tau_1_Sm = 1.5*(t_2_Sm-t_1_Sm);
t_d_Sm = t_2_Sm-tau_1_Sm;
den_Sm = [tau_1_Sm 1];
H_Sm = tf(num, den_Sm)*exp(-t_d_Sm*s);

%f - Metodo de Haglund
t_2_Hag = 15;
tau_1_Hag = t_2_Hag-t_d_approx;
t_d_Hag = t_d_approx;
den_Hag = [tau_1_Hag 1];
H_Hag = tf(num, den_Hag)*exp(-t_d_Hag*s);

%g - Metodo de Sundaresan e Krishnaswami;
y_353_approx = 2 + 0.353*(8-2);
y_853_approx = 2 + 0.853*(8-2);
t_1_Sun = 8.9;
t_2_Sun = 23.3;
tau_1_Sun = 0.67*(t_2_Sun-t_1_Sun);
t_d_Sun = 1.3*t_1_Sun-0.29*t_2_Sun;
den_Sun = [tau_1_Sun 1];
H_Sun = tf(num, den_Sun)*exp(-t_d_Sun*s);

%h - Metodo da inclinacao inicial.
tau_1_inc_inicial = 17 - t_d_approx;
t_d_inc_inicial = t_d_approx;
den_inc = [tau_1_inc_inicial 1];
H_inc = tf(num, den_inc)*exp(-t_d_inc_inicial*s);

%
tau_1_approx = [tau_1_Nish,tau_1_1const_temp,tau_1_4const_temp,tau_1_ZN,tau_1_Sm,tau_1_Hag,tau_1_Sun,tau_1_inc_inicial]'
t_d_approx_all = [t_d_Nish,t_d_1const_temp,t_d_4const_temp,t_d_ZN,t_d_Sm,t_d_Hag,t_d_Sun,t_d_inc_inicial]'
%%
T_Nish = 2+lsim(H_Nish,v,t);
T_1const = 2+lsim(H_1const,v,t);
T_4const = 2+lsim(H_4const,v,t);
T_ZN = 2+lsim(H_ZN,v,t);
T_Sm = 2+lsim(H_Sm,v,t);
T_Hag = 2+lsim(H_Hag,v,t);
T_Sun = 2+lsim(H_Sun,v,t);
T_inc = 2+lsim(H_inc,v,t);

figure('color',[1 1 1])
plot(t,T)
hold all
grid on
plot(t,T_Nish)
plot(t,T_1const)
plot(t,T_4const)
plot(t,T_ZN)
plot(t,T_Sm)
plot(t,T_Hag)
plot(t,T_Sun)
plot(t,T_inc)
legend('Real system','Nish','1 const.','4 const.','Ziegler-Nichols','Smith','Haglund','Sundaresan','Inc')
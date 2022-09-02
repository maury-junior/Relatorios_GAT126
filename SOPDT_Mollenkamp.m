function H = SOPDT_Mollenkamp(K,t1,t2,t3)

% t1 - 0.15 y(inf)
% t2 - 0.45 y(inf)
% t3 - 0.75 y(inf)
s = tf('s');
x = (t2-t1)/(t3-t1);
zeta = (0.0805-5.547*(0.475-x)^2)/(x-0.356);
if zeta <1
   f1=0.708*2.811^zeta ;
elseif zeta >=1
    f1=2.6*zeta-0.6;
end

wn = f1/(t3-t1);

f2=0.922*1.66^zeta;
theta = t2-f2/wn;

if zeta>=1
    tau_1 = (zeta+sqrt(zeta^2-1))/wn;
    tau_2 = (zeta-sqrt(zeta^2-1))/wn;
end

H = tf(K*wn^2,[1 2*zeta*wn wn^2])*exp(-theta*s);
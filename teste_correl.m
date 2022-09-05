N=512;
u=prbs(N,11,1);
u=u-0.5;
lu=length(u);
y=dlsim([1 0.5],[1 -1.5 0.7],u);
e=randn(N,1);
ye = y+e;

for k=1:length(u)-1
    ruu(k) = autoCorrel(u,k);
    ruy(k) = crossCorrel(u,ye,k);
end
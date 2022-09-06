N=512;
u=prbs(N,11,1);
u=u-0.5;
lu=length(u);
y=dlsim([1 0.5],[1 -1.5 0.7],u);
e=randn(N,1);
ye = y+e;

for k=1:N/2
    ruu_v(k) = autoCorrel(u,k);
    ruy(k) = crossCorrel(u,ye,k);
end

ruu = ruu_v';
for k=1:N/2-1
    ruu = [ruu [zeros(k,1); ruu_v(1:length(ruu)-k)']];
end

H = inv(ruu)*ruy';

figure
plot(H)
dados = load('saida_prbs.txt');
t=dados(:,1);
u=dados(:,2);
y=dados(:,3);

Y = y(1:4:end);
U = u(1:4:end);

phi = [Y(1:end-1) U(1:end-1)];

teta = pinv(phi)*Y(2:end);
A = teta(1);
B = teta(2);
Ya(1) = Y(1);
for k=2:length(Y)
    Ya(k) = A*Y(k-1) + B*U(k-1);
end

figure
plot(Y,'b')
hold on 
plot(Ya,'r')
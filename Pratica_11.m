load('ensaio_prbs.txt')

t = ensaio_prbs(:,1);
u = ensaio_prbs(:,2);
y = ensaio_prbs(:,3);

figure('color',[1 1 1])
plot(t,u)
dados = load('dadosordem2aleatorio.txt');
t = dados(:,1);
v = dados(:,2);
T = dados(:,3);
lv = length(v);
% 
% figure
% subplot(211)
% plot(t,T)
% subplot(212)
% plot(t,v)

U=v;
for i=1:lv-1
    U=[U [zeros(i,1); v(1:lv-i)]];
end

h=inv(U)*T;

plot(t,h)
grid on
axis([0 5 -10 10])

u2 = ones(lv,1);
U2 = u2;
for i=1:lv-1
    U2=[U2 [zeros(i,1); u2(1:lv-i)]];
end
y2 = U2*h;
figure
plot(t,y2)
axis([0 15 0 15])
%%
U=v(1:100);
for i=1:99
    U=[U [zeros(i,1); v(1:100-i)]];
end

h=inv(U)*T(1:100);

plot(t(1:100),h)
axis([0 5 -10 10])
%%
u2 = ones(150,1);
U2 = u2;
for i=1:149
    U2=[U2 [zeros(i,1); u2(1:150-i)]];
end
y2 = U2*h;
figure
plot(t(1:150),y2)
axis([0 8 0 15])
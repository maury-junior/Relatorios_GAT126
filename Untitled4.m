load('PULSO10.DAT')
t = PULSO10(:,1);
vazao = PULSO10(:,2);
nivel = PULSO10(:,3);

% filtragem
for i=1:50
   v(i)=mean(vazao(i:i+8));
end;
for i=51:60
   v(i)=vazao(i);
end;
for i=61:length(vazao)-8
   v(i)=mean(vazao(i:i+8));
end;

% sinal de nivel
% Na realidade este deve ser o sinal que foi para o inversor
% a divisao por 250 muda a escala de 1 a 5V para 4 a 20mA.
n=nivel(2:731)*1000/250;
plot(t(2:731),vazao(2:731),t(2:731),v(1:730));
subplot(211)
plot(t(1:730),n(1:730),'b');
axis([1 730 14 21]);
ylabel('mA');
title('(a)')
subplot(212)
plot(t(2:731),v(1:730),'b');
ylabel('volts');
xlabel('t(s)');
title('(b)')
axis([1 730 3.6 4.0]);

h=1.0440; % tempo de amostragem
d=1; % fator de dizimacao
fim=2^14-730+49;
wv=[v(50:730) ones(1,fim)*3.83];
wn=[n(47:727)' ones(1,fim)*3.78*1000/250];

V=fft(wv);
N=fft(wn);
t=0:d*1.0440:1.0440*(length(wv)-1);
f = 0:1/t(length(t)):1/(d*h);
f=2*pi*f;

figure(1);
subplot(211);
semilogx(f,abs(N),'b-');
title('(a)');
ylabel('|U(jw)|');
subplot(212);
semilogx(f,abs(V),'b-');
title('(b)');
ylabel('|Y(jw)|');
xlabel('w (rad/s)');
load pulso10.dat;
t=pulso10(:,1);
vazao=pulso10(:,2);
nivel=pulso10(:,3);

fmax=0.015;
ind=find(f>=fmax);
xind=ind(1);
fft_v=fft(vazao)
V(1:xind);
fft_n=N(1:xind);
f=f(1:xind);

fft_tf = fft(vazao)./fft(nivel*1000/250);

d=1;
h = 1.0440;
t=0:d*1.0440:1.0440*(length(wv)-1);
f = 0:1/t(length(t)):1/(d*h);
f=2*pi*f;

fft_tfa = abs(fft_tf);
fft_tff = angle(fft_tf);

figure(1);
subplot(211);
semilogx(f, 20*log10(fft_tfa), 'b');
ylabel('|H(jw)| (dB)');
title('(a)');
subplot(212);
semilogx(f, (180/pi)*fft_tff, 'b');
xlabel('w (rad/s)');
title('(b)');
ylabel('fase de H(jw) (graus)');
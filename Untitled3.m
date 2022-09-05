dadosordem2rand = load('dadosordem2aleatorio.txt');
tempo = dadosordem2rand(:,1);
entrada = dadosordem2rand(:,2);
saida = dadosordem2rand(:,3);

plot(tempo,saida)
grid on
hold on

%Em relação a resposta em frequência, eu fiz:
for k=1:2000
    if k>1
        entrada_d(k)=entrada(k-1)-entrada(k);
        saida_d(k)=saida(k-1)-saida(k);
    end
end

H = fft(saida)./fft(entrada);
H_d= fft(saida_d)./fft(entrada_d);

freq = 1/length(tempo)*(0:length(tempo)/2);

figure
subplot(211)
semilogx(2*pi*freq(1:length(freq)),20*log10(abs(H(1:length(freq)))),'k--');
title('(a)');
ylabel('|H(j\omega)| (dB)');
subplot(212)
semilogx(2*pi*freq,angle(H(1:length(freq)))*180/pi,'k--');
title('(b)');
ylabel('fase de H(j\omega) (graus)');
xlabel('\omega (rad/s)')

figure
subplot(211)
semilogx(2*pi*freq(1:length(freq)),20*log10(abs(H_d(1:length(freq)))),'k--');
title('(a)');
ylabel('|H(j\omega)| (dB)');
subplot(212)
semilogx(2*pi*freq,angle(H_d(1:length(freq)))*180/pi,'k--');
title('(b)');
ylabel('fase de H(j\omega) (graus)');
xlabel('\omega (rad/s)')



%%
v_d_fft = fft(v_rand_d);
T_d_fft = fft(T_rand_d);
H = T_d_fft./v_d_fft;

% vetor de frequencias
freq=1/(length(T_rand))*(0:length(T_rand)/2);

figure
subplot(211)
semilogx(2*pi*freq,20*log10(abs(H(1:length(freq)))),'k--');
title('(a)');
ylabel('|H(j\omega)| (dB)');
subplot(212)
semilogx(2*pi*freq,angle(H(1:length(freq)))*180/pi,'k--');
title('(b)');
ylabel('fase de H(j\omega) (graus)');
xlabel('\omega (rad/s)')
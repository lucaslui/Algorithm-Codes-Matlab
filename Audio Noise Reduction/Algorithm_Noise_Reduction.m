function [] = Algorithm_Noise_Reduction(arquivo_sinal, arquivo_noise, gerar_arquivo_saida)

% Script que recebe o arquivo .wav ruidoso, toca o som do mesmo, aplica
% o algoritmo de subtração espectral, cria um .wav com a saída do algoritmo
% e cria um gráfico dos espectros dos sinais.

%==========================================================================
%==========================================================================
% Função wavread() lê o sinal de áudio .wav e retorna um vetor com cada
% valor do sinal e a frequência de amostragem junto com o númuero de bits,
% na qual foi codificada.

[sinal,fs,nbits] = wavread(arquivo_sinal); 
[noise,fsn,nbitsn] = wavread(arquivo_noise); 

%==========================================================================
%==========================================================================
% Rotina que aplica os algoritmos de redução de ruído, no caso:
% - A função SSBoll79() aplica a Subtração Espectral proposto por Boll.
% - A função WienerScalart96() aplica o Filtro de Wiener proposto por
%   Scalart.
% - A função MMSESTSA84() aplica a técnica MMSE baseado em STSA.
% Obs: Scripts encontrados na comunidade do Matlab.

output_SS = SSBoll79(sinal,fs); 
output_Wiener = WienerScalart96(sinal,fs);
output_MMSE = MMSESTSA84(sinal,fs);

%==========================================================================
%==========================================================================
% Trecho da função responsável por pegar o sinal filtrado e criar um 
% arquivo .wav com o mesmo nome do arquivo original mas com o adento da 
% descrição de filtrado. Função wavwrite() cria o arquivo .wav, função 
% strcar() concate strings, função fileparts() retorna dados do arquivo, 
% como extensão, nome e caminho.

if (strcmp(gerar_arquivo_saida,'sim') == 1)
    [pathstr,name,ext] = fileparts(arquivo_sinal);
    wavwrite(output_SS,fs,16,strcat(name,' - Filtrado por SS.wav'));
    wavwrite(output_Wiener,fs,16,strcat(name,' - Filtrado por Wiener.wav'));
    wavwrite(output_MMSE,fs,16,strcat(name,' - Filtrado por MMSE.wav'));
end

%==========================================================================
%==========================================================================
% Função subplot() e stem plota os graficos de ambos os sinais na mesma 
% com o intervalo de valores desejado pela variável t.

N = length(sinal);
M = length(output_SS);
K = length(output_Wiener);
J = length(output_MMSE);

RN = 0 : N-1;
RM = 0 : M-1;
RK = 0 : K-1;
RJ = 0 : J-1;

RN_TP = RN*(1/fs);
RM_TP = RM*(1/fs);
RK_TP = RK*(1/fs);
RJ_TP = RJ*(1/fs);

fig_1 = figure(1);
set(fig_1,'Units', 'normalized', 'Position', [0,0,1,1]);

subplot(2,4,1); plot(RN_TP(1:N-1),sinal(1:N-1),'k'); 
xlabel('Tempo (s)'); 
ylabel('Magnitude'); 
title('(a)'); 
axis tight

subplot(2,4,2); plot(RM_TP(1:M-1),output_SS(1:M-1),'b'); 
xlabel('Tempo (s)'); 
ylabel('Magnitude'); 
title('(b)'); 
axis tight

subplot(2,4,3); plot(RK_TP(1:K-1),output_Wiener(1:K-1),'g'); 
xlabel('Tempo (s)'); 
ylabel('Magnitude'); 
title('(c)'); 
axis tight

subplot(2,4,4); plot(RJ_TP(1:J-1),output_MMSE(1:J-1),'r'); 
xlabel('Tempo (s)'); 
ylabel('Magnitude'); 
title('(d)'); 
axis tight

fft_sinal = abs(fft(sinal));
fft_output_SS = abs(fft(output_SS));
fft_output_Wiener = abs(fft(output_Wiener));
fft_output_MMSE = abs(fft(output_MMSE));

N_2 = ceil(N/2);
M_2 = ceil(M/2);
K_2 = ceil(K/2);
J_2 = ceil(J/2);

RN_Hz = RN*fs/N;
RM_Hz = RM*fs/M;
RK_Hz = RK*fs/K;
RJ_Hz = RJ*fs/J;

subplot(2,4,5); plot(RN_Hz(1:N_2), fft_sinal(1:N_2),'k'); 
xlabel('Frequência (Hz)'); 
ylabel('Magnitude'); 
title('(e)'); 
axis tight

subplot(2,4,6); plot(RM_Hz(1:M_2), fft_output_SS(1:M_2),'b'); 
xlabel('Frequência (Hz)'); 
ylabel('Magnitude'); 
title('(f)'); 
axis tight

subplot(2,4,7); plot(RK_Hz(1:K_2), fft_output_Wiener(1:K_2),'g'); 
xlabel('Frequência (Hz)'); 
ylabel('Magnitude'); 
title('(g)'); 
axis tight

subplot(2,4,8); plot(RJ_Hz(1:J_2), fft_output_MMSE(1:J_2),'r'); 
xlabel('Frequência (Hz)'); 
ylabel('Magnitude'); 
title('(h)'); 
axis tight

%==========================================================================
%==========================================================================
% Trecho da função responsável por criar o espectograma dos sinais.

window_spec = hamming(512);
noverlap = 256;
nfft = 1024;

fig_2 = figure(2);
set(fig_2,'Units', 'normalized', 'Position', [0,0,1,1]);

[S,F,T,P] = spectrogram(sinal, window_spec, noverlap, nfft, fs, 'yaxis');
subplot(2,2,1); 
surf(T,F,10*log10(P),'edgecolor','none'); 
title('(a)'); 
axis tight; 
view(0,90);
set(gca,'clim',[-80 -30]);
xlabel('Tempo (s)'); 
ylabel('Frequência (Hz)');
c = colorbar;
ylabel(c,'Potência');

[S_SS,F_SS,T_SS,P_SS] = spectrogram(output_SS, window_spec, noverlap, nfft, fs, 'yaxis');
subplot(2,2,2); 
surf(T_SS,F_SS,10*log10(P_SS),'edgecolor','none');
title('(b)'); 
axis tight;
view(0,90);
set(gca,'clim',[-80 -30]);
xlabel('Tempo (s)'); 
ylabel('Frequência (Hz)');
c = colorbar;
ylabel(c,'Potência');

[S_WI,F_WI,T_WI,P_WI] = spectrogram(output_Wiener, window_spec, noverlap, nfft, fs, 'yaxis');
subplot(2,2,3);
surf(T_WI,F_WI,10*log10(P_WI),'edgecolor','none');
title('(c)'); 
axis tight; 
view(0,90);
set(gca,'clim',[-80 -30]);
xlabel('Tempo (s)'); ylabel('Frequência (Hz)');
c = colorbar;
ylabel(c,'Potência');

[S_MM,F_MM,T_MM,P_MM] = spectrogram(output_MMSE, window_spec, noverlap, nfft, fs, 'yaxis');
subplot(2,2,4);
surf(T_MM,F_MM,10*log10(P_MM),'edgecolor','none');
title('(d)'); 
axis tight; 
view(0,90);
set(gca,'clim',[-80 -30]);
xlabel('Tempo (s)'); ylabel('Frequência (Hz)');
c = colorbar;
ylabel(c,'Potência');

colormap(hot); 

%==========================================================================
%==========================================================================
% Rotina da função responsável por achar a SSNR (Relação Sinal-Ruído
% Segmentada.

snr_SS = snr_Lui_2(noise(1:M), output_SS);
snr_Wiener = snr_Lui_2(noise(1:K), output_Wiener);
snr_MMSE = snr_Lui_2(noise(1:J), output_MMSE); 

ssnr_SS = segsnr_Lui(noise(1:M), output_SS, fs);
ssnr_Wiener = segsnr_Lui(noise(1:K), output_Wiener, fs);
ssnr_MMSE = segsnr_Lui(noise(1:J), output_MMSE, fs);
 
fprintf('SNR(SS): %0.2f dB\n', snr_SS);
fprintf('SNR(Wiener): %0.2f dB\n', snr_Wiener);
fprintf('SNR(MMSE): %0.2f dB\n', snr_MMSE);

fprintf('SSNR(SS): %0.2f dB\n', ssnr_SS);
fprintf('SSNR(Wiener): %0.2f dB\n', ssnr_Wiener);
fprintf('SSNR(MMSE): %0.2f dB\n', ssnr_MMSE);

%==========================================================================
%==========================================================================

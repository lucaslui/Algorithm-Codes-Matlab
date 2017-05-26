%==========================================================================
% Aluno: Lucas Lui Motta      RA: 192701
% Data: 02/05/2017
% Disciplina: Introd. Comunicações Digitais I
% Pós-graduação: Mestrado
%==========================================================================
% Programa MATLAB para 4-QAM:
% 
% Tx:
%       - Gerar 4-QAM;
%       - Gerar g(t) cosseno levantado;
%       - Transladar o sinal para RF;
% 
% Canal:
%       - Inserir Ruído AWGN;
% Rx:
%       - Transladar para banda base e filtrar usando filtro casado;
%       - Obter o sinal em fase e quadratura;
%       - Amostragem;
%       - Obter diagrama de constelação recebido;
%       - Estimar os símbolos transmitidos;
%
% Calcular a SER (Symbol Error Rate).
%==========================================================================

clc;
clear all;
close all;

%==========================================================================
% Parâmetros gerais do Tx-Rx 4-QAM

% Parâmetros Ajustáveis (F_bits e K com problema):
F_bits = 20; % Criado um fator de bits que serão transmitidos (p/ truncamento).
K = 10; % Fator de super-amostragem para a convolução do sinal com o filtro formatador.
Fs = 1e3; % Bit Rate [bits/segundo] ou Frequência de Amostragem - Ex: 1.000 [Hz]
Fc = 1e6; % Frequência da portadora - Ex: 1 [Mhz].
amplt = 2; % Amplitude do simbolo modulado, que causa impacto direto na energia do sinal.
alfa = 0.8; % Fator alfa cosseno levantado utilizado no filtro formatador.
SNR = 5; % Relação sinal ruído, entre o sinal RF que será transmitido e o ruído AWGN recebido.

% Parâmetros Automáticos:
M_Order = 4; % Ordem da modulação QAM - Ex: 4-QAM.
B_symb = log2(M_Order); % Número de bits por simbolo - Ex: log2(4) = 2 [bits/simbolo].
N_symb = F_bits*B_symb; % Número de simbolos que serão transmitidos - Ex: 10*(2) = 20 [simbolos].
N_bits = B_symb*N_symb; % Número de bits que serão enviados.
Ts = 1/Fs; % Período entre bits chegando na stream.
Fb = Fs/M_Order; % Baud Rate [símbolos/segundo] - Ex: 1.000/4-QAM = 250Hz [simbolos/segundo].
Tb = 1/Fb; % Período de símbolos transmitidos.

%==========================================================================
%                           Data - Stream
% Gerador pseudo-randomico de dados binários para simular a stream a ser
% transmitida com a determinada Fs escolhida.

data = randi([0,1], 1, N_bits);

t = linspace(0+(1/Fs), N_bits/Fs, N_bits);

subplot(4,2,[1,2]); stem(t, data);
title('Stream a ser Transmitida (Gerador Pseudo-Randomico)'); 
axis([0 N_bits/Fs 0 1.3]);

%==========================================================================
%                              Slicer
% Dividindo os dados da stream binária de entrada em dois caminhos,
% aqueles que vão na parte InPhase do sinal e aquelas que vão em 
% Quadrature.

data_I = data(1:2:end);
data_Q = data(2:2:end);

t_I = linspace(0+(1/Fs), (N_bits/2)/Fs, N_bits/2);
t_Q = linspace(0+(1/Fs), (N_bits/2)/Fs, N_bits/2);

subplot(4,2,3); stem(t_I, data_I, 'g');
title('Stream InPhase'); 
axis([0 (N_bits/2)/Fs 0 1.3]);

subplot(4,2,4); stem(t_Q, data_Q, 'r'); 
title('Stream Quadrature'); 
axis([0 (N_bits/2)/Fs 0 1.3]);

%==========================================================================
%                          Mapper-Encoder
% Mapeando os valores da stream de dados de cada parte em valores reais, o                  
% 4-QAM (Quadrature Amplitude Modulation) equivale a 2-PAM (Pulse Amplitude 
% Modulation), um InPhase e outro em Quadrature.

mapper = [-amplt  amplt]; Ymax = amplt; Ymin = -amplt;

data_mapped_I = mapper(1 + data_I(1:end));
data_mapped_Q = mapper(1 + data_Q(1:end));

t_mapI = linspace(0+(1/Fs), (N_bits/2)/Fs, length(data_mapped_I));
t_mapQ = linspace(0+(1/Fs), (N_bits/2)/Fs, length(data_mapped_Q));

subplot(4,2,5); stem(t_mapI, data_mapped_I, 'g');
title('Sinal PAM InPhase'); 
axis([0 (N_bits/2)/Fs Ymin-0.3 Ymax+0.3]);

subplot(4,2,6); stem(t_mapQ, data_mapped_Q, 'r');
title('Sinal PAM Quadrature');
axis([0 (N_bits/2)/Fs Ymin-0.3 Ymax+0.3]);

%==========================================================================
%                       Pulse Shaping Filter 
% Após o mapeamento da stream em simbolos a serem transmitidos, os símbolos 
% são passados por um filtro formatador. O sinal utilizado como pulso
% formatador é o frequentemente utilizado cosseno levantado.

rcos_coefs = rcosfir(alfa, 3, K, Tb, 'sqrt');

signal_data_I = conv(rcos_coefs, upsample(data_mapped_I, K));
signal_data_Q = conv(rcos_coefs, upsample(data_mapped_Q, K));

signal_data_I = signal_data_I(31:end-30); % o filtro gera uma cauda de tamanho:
signal_data_Q = signal_data_Q(31:end-30); % ((2*3)+1)*K = 70, que é removida.

t_filt = linspace(0+(1/Fs), (N_bits/2)/Fs, length(signal_data_Q));

figure(1);
subplot(4,2,7); plot(t_filt, signal_data_I, 'g');
title('Sinal InPhase contínuo pós-filtro'); 
axis([0 (N_bits/2)/Fs Ymin Ymax]);

subplot(4,2,8); plot(t_filt, signal_data_Q', 'r');
title('Sinal Quadrature contínuo pós-filtro'); 
axis([0 (N_bits/2)/Fs Ymin Ymax]);

%==========================================================================
%                           Up-Converter 
% Move-se o sinal de banda base para banda passante, e soma ambos sinais 
% InPhase e Quadrature para transmissão pela antena.

t_uc = linspace(0, (N_bits/2)/Fs, length(signal_data_I));

Tx_I =  signal_data_I.*(sqrt(2)*cos(2*pi*Fc*t_uc));
Tx_Q =  signal_data_Q.*(sqrt(2)*sin(2*pi*Fc*t_uc));

signal_RF = Tx_I - Tx_Q;

t_rf = linspace(0+(1/Fs), (N_bits/2)/Fs, length(signal_RF));

figure(2);
subplot(4,1,1); plot(t_rf, Tx_I, 'g');
title('Sinal em Banda Passante (InPhase)'); axis tight;

subplot(4,1,2); plot(t_rf, Tx_Q, 'r');
title('Sinal em Banda Passante (Quadrature)'); axis tight;

subplot(4,1,3); plot(t_rf, Tx_I, 'g', t_rf, Tx_Q, 'r');
title('Ambos sinais Banda Passante (InPhase e Quadrature)'); axis tight;

subplot(4,1,4); plot(t_rf, signal_RF, 'k');
title('Sinal RF - Antena'); axis tight;

%==========================================================================
%                           Signal Corrupted
% Na transmissão do sinal no canal, o mesmo é prejudicado pela presença de
% ruído, para simular esse ruído aleatório é adotado um ruído do tipo
% branco gaussiano (AWGN).

signal_RFandNoise = awgn(signal_RF, SNR);
%signal_RFandNoise = signal_RF;

figure(3);
subplot(6,1,1); plot(t_rf, signal_RF, 'k');
title('Sinal RF - Transmitido'); axis tight;

subplot(6,1,2); plot(t_rf, signal_RFandNoise, 'm');
title('Sinal RF + Ruído AWGN c/ SNR 5 - Recebido'); axis tight;

%==========================================================================
%                           Down-Converter
% O primeiro passo do receptor é mover o sinal na banda passante para banda
% base utilizando o down-converter, que é formado pelo filtro U(f) e o
% o mixer com a função de deslocamento.

baseBand_Rx_I =  signal_RFandNoise.*sqrt(2).*cos(2*pi*Fc*t_uc);
baseBand_Rx_Q =  signal_RFandNoise.*-sqrt(2).*sin(2*pi*Fc*t_uc);

Fpass = 3000; Fstop = 3750; Apass = 1; Astop = 40; Fsamp = 10000;    
obj_h = fdesign.lowpass('fp,fst,ap,ast', Fpass, Fstop, Apass, Astop, Fsamp);    
Hd = design(obj_h, 'equiripple');

baseBand_Rx_I_int = [baseBand_Rx_I zeros(1,10)];
baseBand_Rx_Q_int = [baseBand_Rx_Q zeros(1,10)];

baseBand_Rx_I_LPF = filter(Hd, baseBand_Rx_I_int);
baseBand_Rx_Q_LPF = filter(Hd, baseBand_Rx_Q_int);

baseBand_Rx_I_LPF = baseBand_Rx_I_LPF(11:end);
baseBand_Rx_Q_LPF = baseBand_Rx_Q_LPF(11:end);

subplot(6,1,3); plot(t_rf, baseBand_Rx_I, 'g');
title('Sinal InPhase deslocado p/ Banda Base'); axis tight;

subplot(6,1,4); plot(t_rf, baseBand_Rx_I_LPF, 'g');
title('Sinal InPhase pós-LPF'); axis tight;

subplot(6,1,5); plot(t_rf, baseBand_Rx_Q, 'r');
title('Sinal Quadrature deslocado p/ Banda Base'); axis tight;

subplot(6,1,6); plot(t_rf, baseBand_Rx_Q_LPF, 'r');
title('Sinal Quadrature pós-LPF'); axis tight;

%==========================================================================
%                         Matched Filter
% Utilizando filtro casado para redução do ruído causado no sinal.

signal_Rx_I = conv(rcos_coefs, baseBand_Rx_I_LPF);
signal_Rx_Q = conv(rcos_coefs, baseBand_Rx_Q_LPF);

signal_Rx_I = signal_Rx_I(31:end-30); 
signal_Rx_Q = signal_Rx_Q(31:end-30);

figure(4);
subplot(4,1,1); plot(t_rf, baseBand_Rx_I_LPF, 'g');
title('Sinal InPhase Ruidoso'); axis tight;

subplot(4,1,2); plot(t_rf, signal_Rx_I, 'g');
title('Sinal InPhase pós-Filtro Casado'); axis tight;

subplot(4,1,3); plot(t_rf, baseBand_Rx_Q_LPF, 'r');
title('Sinal Quadrature Ruidoso'); axis tight;

subplot(4,1,4); plot(t_rf, signal_Rx_Q, 'r');
title('Sinal Quadrature pós-Filtro Casado'); axis tight;

%==========================================================================
%                             Sampling
% Amostrando o sinal

data_Rx_I_sampl = signal_Rx_I(1:length(signal_Rx_I)/N_symb:end);
data_Rx_Q_sampl = signal_Rx_Q(1:length(signal_Rx_Q)/N_symb:end);

t_I = linspace(0+(1/Fs), (N_bits/2)/Fs, N_symb);
t_Q = linspace(0+(1/Fs), (N_bits/2)/Fs, N_symb);

figure(5);
subplot(4,2,1); plot(t_rf, signal_Rx_I, 'g');
title('Sinal InPhase pós-Filtro Casado'); axis tight;

subplot(4,2,2); plot(t_rf, signal_Rx_Q, 'r');
title('Sinal Quadrature pós-Filtro Casado'); axis tight;

subplot(4,2,3); stem(t_I, data_Rx_I_sampl, 'g');
title('Sinal InPhase Amostrado'); 
axis([0 (N_bits/2)/Fs Ymin-1 Ymax+1]);

subplot(4,2,4); stem(t_Q, data_Rx_Q_sampl, 'r'); 
title('Sinal Quadrature Amostrado'); 
axis([0 (N_bits/2)/Fs Ymin-1 Ymax+1]);

%==========================================================================
%                  Comparator - "Demapper" - Decoding
% Essa operação trunca o valor recebido no simbolo mais próximo
% definido anteriormente no mapper do transmissor.

data_Rx_I = zeros(1, N_bits/2);
data_Rx_Q = zeros(1, N_bits/2);

for i = 1:length(data_Rx_I_sampl)
    if data_Rx_I_sampl(i) <= 0
        data_Rx_I(i) = 0;
    else
        data_Rx_I(i) = 1;
    end
end

for i = 1:length(data_Rx_Q_sampl)
    if data_Rx_Q_sampl(i) <= 0
        data_Rx_Q(i) = 0;
    else
        data_Rx_Q(i) = 1;
    end
end

subplot(4,2,5); stem(t_I, data_Rx_I, 'g');
title('Stream InPhase - Descodificado'); 
axis([0 (N_bits/2)/Fs 0 1.3]);

subplot(4,2,6); stem(t_Q, data_Rx_Q, 'r'); 
title('Stream Quadrature - Descodificado'); 
axis([0 (N_bits/2)/Fs 0 1.3]);

%==========================================================================
%                Grouping signals InPhase and Quadrature
% Formando a strem de dados final, agrupando 1 a 1 os dados que estavam 
% InPhase e em QuadratureAmostrando:

data_Rx = zeros(N_bits,1);

data_Rx(1:2:end) = data_Rx_I;
data_Rx(2:2:end) = data_Rx_Q;

t = linspace(0+(1/Fs), N_bits/Fs, N_bits);

subplot(4,2,[7,8]); stem(t, data_Rx, 'b');
title('Stream Recebida'); 
axis([0 N_bits/Fs 0 1.3]);

%==========================================================================
%                          Spectro - (Optional) 
% Espectro dos sinais na banda base e banda passante.

spectro_Tx_I = abs(fftshift(fft(signal_data_I))/length(signal_data_I));
spectro_Tx_Q = abs(fftshift(fft(signal_data_Q))/length(signal_data_Q));

spectro_RF = abs(fft(signal_RF)/length(signal_RF));

spectro_Rx_I_LPF = abs(fftshift(fft(baseBand_Rx_I_LPF))/length(baseBand_Rx_I_LPF));
spectro_Rx_Q_LPF = abs(fftshift(fft(baseBand_Rx_Q_LPF))/length(baseBand_Rx_Q_LPF));

spectro_Rx_I_MatchedFilter = abs(fftshift(fft(signal_Rx_I))/length(signal_Rx_I));
spectro_Rx_Q_MatchedFilter = abs(fftshift(fft(signal_Rx_Q))/length(signal_Rx_Q));

f_base = (-length(signal_data_I)/2:(length(signal_data_I)/2)-1)*Fs*K/length(signal_data_I);
%f_pass = (0:(length(spectro_RF))-1)*Fc/length(spectro_RF);
f_pass = linspace(0,Fc*Fs,length(spectro_RF));

figure(6);
subplot(3,3,1); plot(f_base, spectro_Tx_I, 'g');
title('Espectro do sinal InPhase no Transmissor'); axis tight;
subplot(3,3,2); plot(f_base, spectro_Rx_I_LPF, 'g');
title('Espectro sinal ruidoso em Banda Base (InPhase) no Receptor'); axis tight;
subplot(3,3,3); plot(f_base, spectro_Rx_I_MatchedFilter, 'g');
title('Espectro sinal após filtro casado no Receptor'); axis tight;

subplot(3,3,4); plot(f_base, spectro_Tx_Q, 'r');
title('Espectro do sinal Quadrature no Transmissor'); axis tight;
subplot(3,3,5); plot(f_base, spectro_Rx_Q_LPF, 'r');
title('Espectro sinal ruidoso em Banda Base (Quadrature) no Receptor'); axis tight;
subplot(3,3,6); plot(f_base, spectro_Rx_Q_MatchedFilter, 'r');
title('Espectro sinal após filtro casado no Receptor'); axis tight;

subplot(3,3,[7,8,9]); plot(f_pass, spectro_RF, 'k');
title('Especto Sinal RF sendo emitido pelo Transmissor'); axis tight;

%==========================================================================
%                       Constellation Points
% Plotando as constelações dos simbolos que serão transmistidos.

figure(7);
plot(data_mapped_I, data_mapped_Q, 'or', data_Rx_I_sampl, data_Rx_Q_sampl, 'ob');
title('Constellation Points - Tx e Rx'); 
axis([Ymin-0.7 Ymax+0.7 Ymin-0.7 Ymax+0.7]);

%==========================================================================
%                       Symbol Error Rate
% Calculando o número de erros de símbolo e a taxa associas.

[Number_Error, Ratio_Error] = symerr(data,data_Rx');

Ratio_Error = Ratio_Error + 0.0001; % Para evitar erro na hora do plot
                                    % quando o valor de erro é nulo.
figure(8);
pmf_X = [Ratio_Error (1-Ratio_Error)];
pobj = pie(pmf_X);
title('Symbol Error Rate');

hp_color = findobj(pobj, 'Type', 'patch');
set(hp_color(1), 'FaceColor', 'r');
set(hp_color(2), 'FaceColor', 'g');

hp_text = findobj(pobj, 'Type', 'text');
percentValues = get(hp_text,'String');
txt = {'Error: ';'Correct: '}; 
combinedtxt = strcat(txt, percentValues);
oldExtents_cell = get(hp_text,'Extent');
oldExtents = cell2mat(oldExtents_cell);
set(hp_text(1), 'String', combinedtxt(1));
set(hp_text(2), 'String', combinedtxt(2));

%==========================================================================
%                        Algorithm to Save Figures

save_obj_h = get(0,'children');
for i=1:length(save_obj_h)
  saveas(save_obj_h(i), ['figure' num2str(i)], 'jpg');
end

%==========================================================================

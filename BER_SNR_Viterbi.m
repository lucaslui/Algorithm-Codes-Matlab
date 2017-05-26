%==========================================================================
% Aluno: Lucas Lui Motta      RA: 192701
% Data: 07/05/2017
% Disciplina: Introd. Comunicações Digitais I
% Pós-graduação: Mestrado
%==========================================================================
% Programa MATLAB para Análise do Desempenho da modulação com ruído AWGN
% ou seja, simular o BER vs SNR de uma modulação 2-PSK e comparar com valor
% teórico.
% 
% p/ conferência: BER = 10^(-5) com SNR = 9,6dB
%
%==========================================================================

clc;
clear all;
close all;

%==========================================================================

SNRdB = 1:0.5:10;      
SNR = 10.^(SNRdB/10);
Ea = 1;

nBits = 1e6;

mk = [1 0.5 0.25 0.125 0.0625];

BER_simulator = zeros(size(SNR));
BER_noise_ch = zeros(size(SNR));
BER_viterbi = zeros(size(SNR));

a = randsrc(1,nBits, [-sqrt(Ea) sqrt(Ea); 1/2 1/2]);

for i = 1:length(SNR)  

    x = zeros(1,nBits);
    x_noise_ch = zeros(1,nBits);
    
    n = sqrt(Ea/(2*SNR(i)))*randn(1,nBits);
    
    z_awgn = a + n;    
    z_ch_awgn = conv(z_awgn, mk);    
    z_viterbi = real(mlseeq(z_ch_awgn,mk,[-sqrt(Ea);sqrt(Ea)],10,'rst'));
    
    for j = 1:nBits  
        if z_awgn(j) > 0 
            x(j) = sqrt(Ea);
        else
            x(j) = -sqrt(Ea);
        end
    end 
    
    for j = 1:nBits  
        if z_ch_awgn(j) > 0 
            x_noise_ch(j) = sqrt(Ea);
        else
            x_noise_ch(j) = -sqrt(Ea);
        end
    end
    [nErrs_1, BER_simulator(i)] = symerr(x(1:nBits),a);
    [nErrs_2, BER_noise_ch(i)] = symerr(x_noise_ch(1:nBits),a);
    [nErrs_3, BER_viterbi(i)] = symerr(z_viterbi(1:nBits),a);
end

%==========================================================================
%                            Análitico

BER_analytical = qfunc(sqrt(2*SNR));

%==========================================================================
%                              Plot

figure(1);
semilogy(SNRdB, BER_analytical, 'b'); hold on;
semilogy(SNRdB, BER_simulator, 'db'); 
semilogy(SNRdB, BER_noise_ch, 'dr');
semilogy(SNRdB, BER_viterbi, 'dg'); hold off;

title('2-PSK and AWGN Simulation');
xlabel('Eb/No (dB)');
ylabel('BER');
legend('Analytical w/ AWGN','Simulated w/ AWGN', 'Noise Channel + AWGN', 'Used Viterbi');

%==========================================================================
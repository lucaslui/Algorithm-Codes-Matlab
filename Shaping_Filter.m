% Criando um interpolador no MATLAB

%==========================================================================
% Equação que define a interpolação dos símbolos:   

am = [3 5 2 1 4 5 -2 -4 -5 -1]; % Símbolos que seram transmitidos

T = 1;  % Período do baud rate em segundos - ex: T = 1ms

fb = 1/T;   % Frequência do baud rate em hertz - ex: fb = 1kHz 

L = length(am); % Cálculo do tamanho do vetor de simbolos - ex: [3 5 2 1 4 5 -2 -4 -5 -1] = 10

m = 1:L;  % Criação do vetor dos indices das amostras - ex: [1 2 3 4 5 6 7 8 9 10] 

t_am = m*T; % Base de tempo para o sinal discreto dos simbolos - ex: [1*0.001 2*0.001 ... 10*0.001]

fs = 1000; % Frequência de amostragem.

t = linspace(0,L*T,fs); % Eixo Analógico - ex: [0 L/100 (L/100)*2 ... L] = [0 0,05 0,10 ... 5]

% s_inter_sinc_1 = sinc(pi*(t - 1*T)/T)*am(1); % Teste

s_inter_sinc_1 = sinc(t - 1*T)*am(1); % Cada símbolo interpolado com uma função sinc
s_inter_sinc_2 = sinc(t - 2*T)*am(2);
s_inter_sinc_3 = sinc(t - 3*T)*am(3);
s_inter_sinc_4 = sinc(t - 4*T)*am(4);
s_inter_sinc_5 = sinc(t - 5*T)*am(5);
s_inter_sinc_6 = sinc(t - 6*T)*am(6);
s_inter_sinc_7 = sinc(t - 7*T)*am(7);
s_inter_sinc_8 = sinc(t - 8*T)*am(8);
s_inter_sinc_9 = sinc(t - 9*T)*am(9);
s_inter_sinc_10 = sinc(t - 10*T)*am(10);

[TEMPO, MT] = ndgrid(t,m*T); % Ferramenta matlab para subtração dos termos enquanto na somatória.

s_inter_sinc = sinc(TEMPO - MT)*(am');   % Desloc. a função Sinc. e multiplica pelo Símbolo, e soma todas funções Sinc. criadas.
s_inter_gate = (heaviside(TEMPO - MT/T)- heaviside(TEMPO - (MT+1)))*(am'); % Função Gate de intervalo T equivale a: u(n-0) - u(n-T). 
s_inter_reta = am(m); % Função Gate de intervalo T equivale a: u(n-0) - u(n-T). 

%==========================================================================
% Saída do script (gráficos no domínio do tempo) 

fig_1 = figure(1);

subplot(2,2,1); 
plot(t_am, am, 'o', t, s_inter_sinc_1, t, s_inter_sinc_2, t, s_inter_sinc_3, t, s_inter_sinc_4, t, s_inter_sinc_5...
                  , t, s_inter_sinc_6, t, s_inter_sinc_7, t, s_inter_sinc_8, t, s_inter_sinc_9, t, s_inter_sinc_10);
xlabel('Tempo (s)'); 
ylabel('Amplitude'); 
title('Símbolos que serão transmitidos e seus Sinc');
legend('Símbolo', 'Sinc() x Símbolo');
axis tight;

subplot(2,2,2); 
plot(t, s_inter_sinc, 'b', t_am, am, 'or'); 
xlabel('Tempo (s)'); 
ylabel('Amplitude'); 
legend('Interpolado', 'Símbolo');
title('Interpolação utilizando Sinc'); 
axis tight;

subplot(2,2,3); 
plot(t_am, s_inter_reta, 'b', t_am, am, 'or'); 
xlabel('Tempo (s)'); 
ylabel('Amplitude');
legend('Interpolado', 'Símbolo');
title('Interporlação utilizando Retas'); 
axis tight;

subplot(2,2,4); 
plot(t, s_inter_gate, 'b', t_am, am, 'or'); 
xlabel('Tempo (s)'); 
ylabel('Amplitude');
legend('Interpolado', 'Símbolo');
title('Interporlação utilizando Gates'); 
axis tight;

%==========================================================================
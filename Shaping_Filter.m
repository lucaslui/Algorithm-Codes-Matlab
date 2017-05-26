% Criando um interpolador no MATLAB

%==========================================================================
% Equa��o que define a interpola��o dos s�mbolos:   

am = [3 5 2 1 4 5 -2 -4 -5 -1]; % S�mbolos que seram transmitidos

T = 1;  % Per�odo do baud rate em segundos - ex: T = 1ms

fb = 1/T;   % Frequ�ncia do baud rate em hertz - ex: fb = 1kHz 

L = length(am); % C�lculo do tamanho do vetor de simbolos - ex: [3 5 2 1 4 5 -2 -4 -5 -1] = 10

m = 1:L;  % Cria��o do vetor dos indices das amostras - ex: [1 2 3 4 5 6 7 8 9 10] 

t_am = m*T; % Base de tempo para o sinal discreto dos simbolos - ex: [1*0.001 2*0.001 ... 10*0.001]

fs = 1000; % Frequ�ncia de amostragem.

t = linspace(0,L*T,fs); % Eixo Anal�gico - ex: [0 L/100 (L/100)*2 ... L] = [0 0,05 0,10 ... 5]

% s_inter_sinc_1 = sinc(pi*(t - 1*T)/T)*am(1); % Teste

s_inter_sinc_1 = sinc(t - 1*T)*am(1); % Cada s�mbolo interpolado com uma fun��o sinc
s_inter_sinc_2 = sinc(t - 2*T)*am(2);
s_inter_sinc_3 = sinc(t - 3*T)*am(3);
s_inter_sinc_4 = sinc(t - 4*T)*am(4);
s_inter_sinc_5 = sinc(t - 5*T)*am(5);
s_inter_sinc_6 = sinc(t - 6*T)*am(6);
s_inter_sinc_7 = sinc(t - 7*T)*am(7);
s_inter_sinc_8 = sinc(t - 8*T)*am(8);
s_inter_sinc_9 = sinc(t - 9*T)*am(9);
s_inter_sinc_10 = sinc(t - 10*T)*am(10);

[TEMPO, MT] = ndgrid(t,m*T); % Ferramenta matlab para subtra��o dos termos enquanto na somat�ria.

s_inter_sinc = sinc(TEMPO - MT)*(am');   % Desloc. a fun��o Sinc. e multiplica pelo S�mbolo, e soma todas fun��es Sinc. criadas.
s_inter_gate = (heaviside(TEMPO - MT/T)- heaviside(TEMPO - (MT+1)))*(am'); % Fun��o Gate de intervalo T equivale a: u(n-0) - u(n-T). 
s_inter_reta = am(m); % Fun��o Gate de intervalo T equivale a: u(n-0) - u(n-T). 

%==========================================================================
% Sa�da do script (gr�ficos no dom�nio do tempo) 

fig_1 = figure(1);

subplot(2,2,1); 
plot(t_am, am, 'o', t, s_inter_sinc_1, t, s_inter_sinc_2, t, s_inter_sinc_3, t, s_inter_sinc_4, t, s_inter_sinc_5...
                  , t, s_inter_sinc_6, t, s_inter_sinc_7, t, s_inter_sinc_8, t, s_inter_sinc_9, t, s_inter_sinc_10);
xlabel('Tempo (s)'); 
ylabel('Amplitude'); 
title('S�mbolos que ser�o transmitidos e seus Sinc');
legend('S�mbolo', 'Sinc() x S�mbolo');
axis tight;

subplot(2,2,2); 
plot(t, s_inter_sinc, 'b', t_am, am, 'or'); 
xlabel('Tempo (s)'); 
ylabel('Amplitude'); 
legend('Interpolado', 'S�mbolo');
title('Interpola��o utilizando Sinc'); 
axis tight;

subplot(2,2,3); 
plot(t_am, s_inter_reta, 'b', t_am, am, 'or'); 
xlabel('Tempo (s)'); 
ylabel('Amplitude');
legend('Interpolado', 'S�mbolo');
title('Interporla��o utilizando Retas'); 
axis tight;

subplot(2,2,4); 
plot(t, s_inter_gate, 'b', t_am, am, 'or'); 
xlabel('Tempo (s)'); 
ylabel('Amplitude');
legend('Interpolado', 'S�mbolo');
title('Interporla��o utilizando Gates'); 
axis tight;

%==========================================================================
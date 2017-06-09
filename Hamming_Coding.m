% =========================================================================
% CAP 7 - Teoria da Informação
% 
% Código de Hamming (Classe código de Bloco)
% 
% Objetivo: Construir um gráfico da probabilidade de erro por prob. do 
% canal BSC (Pb(alfa)).
% =========================================================================

clc;
clear all;
close all;

% =========================================================================
% =========================================================================
%
%                       Hamming com paridade 2
%
% =========================================================================
% =========================================================================

k_compParidade = 2; % Com bloco do codificador de Hamming(N,n) = Hamm(7,4)
N_compTotal = 2^k_compParidade - 1; % Obs.: (k_compParidade >= 2)
n_compDados = 2^k_compParidade - k_compParidade - 1; 

alfa = 0:0.1:1;

nTest = 1000;
probErro_k2 = zeros(1,length(alfa));

for i = 1:length(alfa)
    
    acumErrs = 0; % Variável para contar o número de erros total por alfa no BSC
    
    for j = 1:nTest
        
        w = randsrc(1,n_compDados,[0 1; 0.5 0.5]); % Dados originais
        
        p_matrix = [1; % Matriz de Paridade
                    1]; 
        h_matrix = [p_matrix eye(N_compTotal - n_compDados)]; % Matriz H = [Paridade I(3X3)]
        g_matrix = [eye(n_compDados) p_matrix']; % Matriz G -> Geradora (TX)
        e_matrix = fliplr(eye(N_compTotal)); % Matriz E -> Indicador do Erro
        s_matrix = mod(e_matrix*h_matrix',2); % Matriz da Síndrome
                
        x = mod(w*g_matrix,2); % Sinal codificado por Hamming que será transmitido
        
        y = bsc(x, alfa(i)); % Sinal recebido após passagem pelo canal BSC
                
        s_signal = mod(y*h_matrix',2); % Detecta se houve erro através do cálculo da síndrome
        
        % Através da sindrome, procura o erro correspondente:
        errorRow = find(ismember(s_matrix, s_signal,'rows'));  
        
        if errorRow ~= 0
            x_estimate = xor(y, e_matrix(errorRow,:)); % Executa a correção
        else
            x_estimate = y; % Caso não haja detecção de erro, não a correção
        end
        
        if biterr(x, x_estimate) ~= 0 % Faz a comparação do sinal estimado no RX com o sinal transmitido.
            acumErrs = acumErrs + 1;  % Caso houve qualquer erro, adiciona a contagem de erros.
        end
    end
    probErro_k2(i) = acumErrs/nTest;
end

% =========================================================================
% =========================================================================
%
%                       Hamming com paridade 3
%
% =========================================================================
% =========================================================================

k_compParidade = 3; % Com bloco do codificador de Hamming(N,n) = Hamm(7,4)
N_compTotal = 2^k_compParidade - 1; % Obs.: (k_compParidade >= 2)
n_compDados = 2^k_compParidade - k_compParidade - 1; 

alfa = 0:0.1:1;

nTest = 100;
probErro_k3 = zeros(1,length(alfa));

for i = 1:length(alfa)
    
    acumErrs = 0; % Variável para contar o número de erros total por alfa no BSC
    
    for j = 1:nTest
        
        w = randsrc(1,n_compDados,[0 1; 0.5 0.5]); % Dados originais
        
        p_matrix = [1 1 0 1; % Matriz de Paridade
                    1 0 1 1;
                    0 1 1 1]; 
        h_matrix = [p_matrix eye(N_compTotal - n_compDados)]; % Matriz H = [Paridade I(3X3)]
        g_matrix = [eye(n_compDados) p_matrix']; % Matriz G -> Geradora (TX)
        e_matrix = fliplr(eye(N_compTotal)); % Matriz E -> Indicador do Erro
        s_matrix = mod(e_matrix*h_matrix',2); % Matriz da Síndrome
                
        x = mod(w*g_matrix,2); % Sinal codificado por Hamming que será transmitido
        
        y = bsc(x, alfa(i)); % Sinal recebido após passagem pelo canal BSC
                
        s_signal = mod(y*h_matrix',2); % Detecta se houve erro através do cálculo da síndrome
        
        % Através da sindrome, procura o erro correspondente:
        errorRow = find(ismember(s_matrix, s_signal,'rows'));  
        
        if errorRow ~= 0
            x_estimate = xor(y, e_matrix(errorRow,:)); % Executa a correção
        else
            x_estimate = y; % Caso não haja detecção de erro, não a correção
        end
        
        if biterr(x, x_estimate) ~= 0 % Faz a comparação do sinal estimado no RX com o sinal transmitido.
            acumErrs = acumErrs + 1;  % Caso houve qualquer erro, adiciona a contagem de erros.
        end
    end
    probErro_k3(i) = acumErrs/nTest;
end

% =========================================================================
% =========================================================================
%
%                       Hamming com paridade 4
%
% =========================================================================
% =========================================================================

k_compParidade = 4; % Com bloco do codificador de Hamming(N,n) = Hamm(7,4)
N_compTotal = 2^k_compParidade - 1; % Obs.: (k_compParidade >= 2)
n_compDados = 2^k_compParidade - k_compParidade - 1; 

probErro_k4 = zeros(1,length(alfa));

for i = 1:length(alfa)
    
    acumErrs = 0; % Variável para contar o número de erros total por alfa no BSC
    
    for j = 1:nTest
        
        w = randsrc(1,n_compDados,[0 1; 0.5 0.5]); % Dados originais
        
        p_matrix = [1 1 0 1 1 0 1 0 1 0 1; % Matriz de Paridade
                    1 0 1 1 0 1 1 0 0 1 1; 
                    0 1 1 1 0 0 0 1 1 1 1; 
                    0 0 0 0 1 1 1 1 1 1 1]; 
        h_matrix = [p_matrix eye(N_compTotal - n_compDados)]; % Matriz H = [Paridade I(3X3)]
        g_matrix = [eye(n_compDados) p_matrix']; % Matriz G -> Geradora (TX)
        e_matrix = fliplr(eye(N_compTotal)); % Matriz E -> Indicador do Erro
        s_matrix = mod(e_matrix*h_matrix',2); % Matriz da Síndrome
                
        x = mod(w*g_matrix,2); % Sinal codificado por Hamming que será transmitido
        
        y = bsc(x, alfa(i)); % Sinal recebido após passagem pelo canal BSC
                
        s_signal = mod(y*h_matrix',2); % Detecta se houve erro através do cálculo da síndrome
        
        % Através da sindrome, procura o erro correspondente:
        errorRow = find(ismember(s_matrix, s_signal,'rows'));  
        
        if errorRow ~= 0
            x_estimate = xor(y, e_matrix(errorRow,:)); % Executa a correção
        else
            x_estimate = y; % Caso não haja detecção de erro, não a correção
        end
        
        if biterr(x, x_estimate) ~= 0 % Faz a comparação do sinal estimado no RX com o sinal transmitido.
            acumErrs = acumErrs + 1;  % Caso houve qualquer erro, adiciona a contagem de erros.
        end
    end
    probErro_k4(i) = acumErrs/nTest;
end

% =========================================================================
% =========================================================================
%
%                       Hamming com paridade 5
%
% =========================================================================
% =========================================================================

k_compParidade = 5; % Com bloco do codificador de Hamming(N,n) = Hamm(7,4)
N_compTotal = 2^k_compParidade - 1; % Obs.: (k_compParidade >= 2)
n_compDados = 2^k_compParidade - k_compParidade - 1; 

probErro_k5 = zeros(1,length(alfa));

for i = 1:length(alfa)
    
    acumErrs = 0; % Variável para contar o número de erros total por alfa no BSC
    
    for j = 1:nTest
        
        w = randsrc(1,n_compDados,[0 1; 0.5 0.5]); % Mensagem 1xtamBloco = [x1 x2 x3 x4]
        
        % Matriz de Paridade
        p_matrix = [1 1 0 1 1 0 1 0 1 0 1 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1; 
                    1 0 1 1 0 1 1 0 0 1 1 0 1 1 0 0 1 1 0 0 1 1 0 0 1 1; 
                    0 1 1 1 0 0 0 1 1 1 1 0 0 0 1 1 1 1 0 0 0 0 1 1 1 1; 
                    0 0 0 0 1 1 1 1 1 1 1 0 0 0 0 0 0 0 1 1 1 1 1 1 1 1
                    0 0 0 0 0 0 0 0 0 0 0 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1]; 
        h_matrix = [p_matrix eye(N_compTotal - n_compDados)]; % Matriz H = [Paridade I(3X3)]
        g_matrix = [eye(n_compDados) p_matrix']; % Matriz G -> Geradora (TX)
        e_matrix = fliplr(eye(N_compTotal)); % Matriz E -> Indicador do Erro
        s_matrix = mod(e_matrix*h_matrix',2); % Matriz da Síndrome
                
        x = mod(w*g_matrix,2); % Sinal codificado por Hamming que será transmitido
        
        y = bsc(x, alfa(i)); % Sinal recebido após passagem pelo canal BSC
                
        s_signal = mod(y*h_matrix',2); % Detecta se houve erro através do cálculo da síndrome
        
        % Através da sindrome, procura o erro correspondente:
        errorRow = find(ismember(s_matrix, s_signal,'rows'));  
        
        if errorRow ~= 0
            x_estimate = xor(y, e_matrix(errorRow,:)); % Executa a correção
        else
            x_estimate = y; % Caso não haja detecção de erro, não a correção
        end
        
        if biterr(x, x_estimate) ~= 0 % Faz a comparação do sinal estimado no RX com o sinal transmitido.
            acumErrs = acumErrs + 1;  % Caso houve qualquer erro, adiciona a contagem de erros.
        end
    end
    probErro_k5(i) = acumErrs/nTest;
end

% =========================================================================

plot(alfa, probErro_k2, 'r', alfa, probErro_k3, 'k', alfa, probErro_k4, 'b', alfa, probErro_k5, 'g');
axis([0 1 0 1.2]);
title('Prob. de Erro c/ Cod. de Canal Hamming');
xlabel('alfa - BSC');
ylabel('Pw - Probabilidade de Erro da Mensagem');
legend('Hamming (3,1)','Hamming(7,4)','Hamming(15,11)', 'Hamming(31,26)');

% =========================================================================
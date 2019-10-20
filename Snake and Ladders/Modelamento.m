% %========================================================================
% % Aluno: Lucas Lui Motta      RA: 192701
% % Data: 05/09/2017
% % Disciplina: Processos Estoc�sticos para Engenharia
% % P�s-gradua��o: Mestrado
% %========================================================================
% % Projeto � Snakes and Ladders
% % 
% % � Plote um gr�fico cujo eixo x seja a rodada e cujo eixo y seja a
% % probabilidade de vencer o jogo naquela rodada.
% % � Plote um gr�fico cujo eixo x seja a rodada e cujo eixo y seja a
% % probabilidade de o jogo acabar ou j� ter acabado naquela rodada.
% % � Prazo: 12 de setembro
% %
% % Dicas: http://www.datagenetics.com/blog/november12011/
% % Fonte: https://www.flickr.com/photos/mythoto/9720925815
% %========================================================================

clc;
clear all;
close all;

%==========================================================================
%                    Simula��o Objetiva (Frequ�ncialista)

num_simulacoes = 100000; % N�mero de simula��es adotado para o experimento.
num_jogadas_max = 1000; % Vari�vel de seguran�a para o experimento.
pmf_win_exp = zeros(1, num_jogadas_max); % Vari�vel que armazenar� a PMF.
cdf_win_exp = zeros(1, num_jogadas_max); % Vari�vel que armazenar� a CDF.

% Definindo casas que cont�m cobras ou escadas. 
escadas = [ 4  6 12 21 25 33 64 67 71 80; 
           38 14 27 78 55 51 94 88 92 99];       
cobras = [ 95 93 83 97 79 77 65 48 41 37;
           74 53 61 29 22 44  2  9 19 15];
% Definindo pontos que ser�o mudados devidos as escadas e cobras.
pulo_start = [escadas(1,:) cobras(1,:)];
pulo_end = [escadas(2,:) cobras(2,:)];

for i = 1:num_simulacoes
    
    venceu = false;  % Vari�vel p/ indicar se o jogador j� venceu.
    estouro = false; % Vari�vel p/ indicar muitas itera��es do experimento.
    num_jogadas = -1; % Valor inicial dado ao n�mero de jogadas, valor -1
                      % pois no la�o interno, existe a jogada 0.
    atual_pos = 0; % Vari�vel que indica a casa atual do jogador.
    
    while (venceu == false) && (estouro == false) 
                
        roll_dado = randi([1 6],1); % Lan�a o dado 1d6 do jogo.
        
        num_jogadas = num_jogadas + 1; % Atualiza n�mero de jogadas 
                                       % at� o momento (inicia em 0).
        
        % Testa se n�o houve estouro.
        if not(num_jogadas >= num_jogadas_max)
            % Atualiza posi��o conforme valor encontrado no dado.
            atual_pos = atual_pos + roll_dado; 
            % La�o de varredura para saber se a casa atual tem alguma 
            % ou cobra, se tiver, atualiza o valor de posi��o atual com o 
            % o valor destino item achado.
            for j = 1:length(pulo_start) 
                if atual_pos == pulo_start(j);
                    atual_pos = pulo_end(j);
                end
            end
            % Executa um teste para saber se o jogador j� venceu ou n�o a
            % partida nessa rodada se sim, sai do la�o.
            if atual_pos >= 100
                venceu = true;
            end 
        else
            estouro = true; % Caso tenha ocorrido estouro, atualiza a var.
        end
    end   
    % Caso n�o tenha ocorrido estouro, � incrementado aquele n�mero de
    % jogadas na vari�vel de contagem, que conta quanto ocorreu de vit�ria
    % em cada valor de rodada.
    if (estouro == false)
        pmf_win_exp(num_jogadas) = pmf_win_exp(num_jogadas) + 1;
    end
end

% Por fim, encontra a probabilidade baseando na divis�o da frequ�ncia de
% ocorr�ncia de vit�ria em cada rodada com  o n�mero total de simula��es.
pmf_win_exp = pmf_win_exp/num_simulacoes;

% Com a PMF calculada acima � poss�vel encontrar a CDF fazendo a soma de 
% todas as probabilidades de cada rodada anterior ou igual i e atribuindo a
% CDF.
for i = 1:length(cdf_win_exp)  
   cdf_win_exp(i) = sum(pmf_win_exp(1:i));
end

%==========================================================================
%                    Simula��o Subjetiva (Bayeseano)

num_casas = 100;    % N�mero de casas do tabuleiro.
num_jogadas = 120;  % N�mero de jogadas analisado.
intervalo = 1:num_jogadas; % Intervalo adotado para analise.

% Inicializando a Matriz de Transi��o (Markov) de 101 linhas x 107 colunas:
% - Com 100 linhas + 1 linha devido a casa 0 (inicial) do tabuleiro.
% - Com 100 colunas + 1 devido vetor matlab e +6 pelas "possiveis casas 
%   sobressalentes" do tabuleiro quando tiver perto do final e o dado
%   1d6 for jogado.
m_trans = zeros(num_casas+1, num_casas+7); 

% Inicializando na Matriz de Transi��o as prob. de transi��o de cada estado
% em rela��o a transi��o da casa i p/ casa j do jogo.
%     j        10    11     12     13    14      15     16     17
%  i [... ... ... ... ... ...... ... ... ... ... ... ... ... ... ...
% 10  ... 0 0.1667 0.1667 0.1667 0.1667 0.1667 0.1667      0      0 ...
% 11  ... 0      0 0.1667 0.1667 0.1667 0.1667 0.1667 0.1667      0 ...
% 12  ... 0      0      0 0.1667 0.1667 0.1667 0.1667 0.1667 0.1667 ...
% 13  ... 0      0      0      0 0.1667 0.1667 0.1667 0.1667 0.1667 ...
%     ... ... ... ... ... ...... ... ... ... ... ... ... ... ...    ...]
for i = 1:(num_casas+1)
    m_trans(i,(i+1:i+6)) = 1/6;
end

% Inicializando na Matriz de Transi��o o valor real das prob. da ultima
% casa do jogo, que cont�m a soma das prob. das casas 101,102,103,104,105.
for i = 1:num_casas+1 
    m_trans(i, num_casas+1) = sum(m_trans(i, num_casas+1:num_casas+7));   
end

% Inicializando na Matriz de Trans. o tamanho correto do n�mero de casas, 
% ou seja, ap�s levado em considera��o as prob. das casas sobressalentes as
% mesmas s�o cortadas da matriz de transi��o p/ 101-101.
m_trans = m_trans(:, 1:num_casas+1); 

% Finalizando a Matriz de Trans. com os valores de cobras e escadas na 
% que nada mais s�o do que altera��es das transi��es dos estados comuns.
%      j     10     11      12      13      14      15     16     
%  i [... ... ... ... ... ...... ... ... ... ... ... ... ... ... ...
% 10  ...      0  0.1667       0  0.1667  0.1667       0  0.1667 ...
% 11  ...      0       0  0.1667  0.1667  0.1667  0.1667  0.1667 ...
% 12  ... 0.1667       0       0  0.1667       0  0.1667  0.1667 ...
% 13  ...      0  0.1667       0       0  0.1667       0  0.1667 ...
%     ... ... ... ... ... ...... ... ... ... ... ... ... ... ... ...]
for i = 1:length(pulo_start)
    altera_casa = m_trans(:, pulo_start(i) + 1);
    indices = find(altera_casa > 0);
    m_trans(indices, pulo_start(i) + 1) = 0;
    m_trans(indices, pulo_end(i) + 1) = m_trans(indices,pulo_end(i) + 1)...
                                        + altera_casa(indices);
end

% Para come�ar a cadeia de Markov � preciso escolher um estado, com as 
% prob. de estar em cada casa. Como nesse jogo se inicia na casa 0 com 
% 100% de certeza de estar na mesma, � criado um vetor linha com 1 na
% posi��o 0 (1 no matlab) e 0 nas demais posi��es.
start_jogo = [1 zeros(1,100)];

% Inicializando com zero os valores que conteram a PMF e a CDF da prob.
% respectivamente vencer o jogo naquela rodada  e do jogo acabar ou j� ter
% acabado naquela rodada. 
pmf_win_sub = zeros(1,num_jogadas);
cdf_win_sub = zeros(1,num_jogadas+1);

% Calculando a prob. acumulada em cada jogada 'i', para isso, � feito a
% atualiza��o das prob. dos estados atrav�s de novas multiplica��es com 
% a Matriz de Transi��o.
for i = 1:length(cdf_win_sub) 
    temp = start_jogo*(m_trans^i);
    cdf_win_sub(i) = temp(num_casas + 1) ;
end

% Com a probabilidade acumulada (CDF) em m�os, pela defini��o da prob.
% acumulada discreta � poss�vel calcular a PMF atrav�s da subtra��o da CDF 
% posterior e a anterior ao instante (no caso, rodada) analisada.
for i = 1:length(pmf_win_sub) 
    pmf_win_sub(i) = cdf_win_sub(i+1) - cdf_win_sub(i);
end

%==========================================================================
%                           Gr�ficos de Sa�da
subplot(2,1,1);
grid on
plot(intervalo, cdf_win_exp(intervalo), 'r'); hold on; grid on;
plot(intervalo, cdf_win_sub(intervalo), 'g'); hold off;
legend('Experimento (Objetivando ou Frequencialista)',... 
       'Modelagem Formal (Subjetivador ou Bayesiano)');
title('Gr�fico da probabilidade de o jogo acabar ou j� ter acabado naquela rodada (Snakes and Ladders)');
xlabel('Rodada'); ylabel('Probabilidade');

subplot(2,1,2);
[max_prob, idx] = max(pmf_win_sub);
h1 = plot(intervalo, pmf_win_exp(intervalo), 'r'); hold on; grid on;
h2 = plot(intervalo, pmf_win_sub(intervalo), 'g'); hold on; 
h3 = plot(intervalo(idx), pmf_win_sub(idx), 'pb'); hold off;
legend('Experimento (Objetivando ou Frequencialista)',... 
       'Modelagem Formal (Subjetivador ou Bayesiano)',...
        sprintf('Maior Probabilidade (Moda) = %0.3f \nRodada = %0.3f', max_prob, idx));
title('Gr�fico da probabilidade de vencer o jogo naquela rodada (Snakes and Ladders)');
xlabel('Rodada'); ylabel('Probabilidade');

% %==========================================================================
% %========================================================================
% % Aluno: Lucas Lui Motta      RA: 192701
% % Data: 16/10/2017
% % Disciplina: Processos Estoc�sticos para Engenharia
% % P�s-gradua��o: Mestrado
% %========================================================================
% %
% % Projeto � Aleatoriedade da Seq. construida por M�quina vs Homem.
% % 
% %========================================================================

clc;
clear all;
close all;

%==========================================================================
%   Simula��o: "The runs test, or the Wald-Wolfowitz test for randomness" 
%
% Baseado: http://www.itl.nist.gov/div898/handbook/eda/section3/eda35d.htm
%          http://blogs.sas.com/content/iml/2013/10/09/how-to-tell-if-a-sequence-is-random.html
%==========================================================================

% Abrindo arquivo de dados com a lista de sequ�ncias geradas.
file = load('BitSequences.mat');
list_seq = file.bits;

% N�vel de signific�ncia utilizado no teste de hipotese.
alfa = 0.1; 

% Encontrando o n�mero de sequ�ncias da lista de sequ�ncias. 
num_seq = length(list_seq(:,1));

cont_alunos = 0; % Contador de sequ�ncias feitas por alunos.
cont_matlab = 0; % Contador de sequ�ncdias feitas pelo Matlab (Aleat�rias).

% Inic. vetor que armazena o resultado do teste de cada sequencia.
vect_result = zeros(1,num_seq);

% Loop que itera sob. todas as seq. realizando o teste de aleatoriedade.
for i = 1:num_seq
    
    % Calculando o n�mero de digitos da sequ�ncia atual.
    tam_seq = length(list_seq(i,:));
    
    % Encontrando a quantidade de 0's e 1's na sequ�ncia.
    n = sum(list_seq(i,:)); % Quantidade de 1's.
    m = (tam_seq - n);      % Quantidade de 0's.

    % Valor esperado, vari�ncia e desvio padr�o do n� de corridas.
    R_med = (2*n*m/(n+m))+1;
    VarR = (2*n*m*((2*n*m)-n-m))/((n+m)^2*(n+m-1));
    SigmaR = sqrt(VarR);

    % Contando o n�mero real de corridas da sequ�ncia atual.
    R = 1;
    for k=2:length(list_seq(i,:))
        if list_seq(i,k)~=list_seq(i,k-1)
            R = R + 1;
        end
    end
    
    % Calculando fator de teste estatistico:
    Z = (R-R_med)/SigmaR;
    
    % Calculando da margem de aceita��o do teste de hip�tese.
    Z_crit = norminv(1-alfa/2)  
  
    % Tomando decis�o baseado nos valores de teste e margem.
    if abs(Z) > Z_crit
        cont_alunos = cont_alunos + 1;
        vect_result(i) = 0;
    else
        cont_matlab = cont_matlab + 1;
        vect_result(i) = 1;        
    end
end

figure(1);
stem(find(vect_result==0), vect_result(find(vect_result==0)), 'r'); hold on;
stem(find(vect_result>0), vect_result(find(vect_result>0)), 'g'); hold on;
legend(sprintf('Pseudo N�o-Aleat�rias (Total - %0.2d)',cont_alunos), ...
       sprintf('Pseudo Aleat�rias (Total - %0.2d)', cont_matlab));
str = [sprintf('Falso Positivo - %0.2d', length(find(vect_result(1:35)==1)))...
       newline ...
       sprintf('Falso Negativo - %0.2d', length(find(vect_result(36:70)==0)))];
annotation('textbox',[.14 .805 .1 .1],'String',str,'FitBoxToText','on');
%title("M�todo de contagem e teste de hip�tese baseado no N�mero de Corridas");
axis([0 71 0 1.5]);


%==========================================================================
%                   Simula��o: "M�todo de Contagem" 
%
% Baseado: http://blog.mrmeyer.com/2011/can-you-recognize-random/

n_bits = 2; % Variavel que define o numero de bits do agrupamento.

cont_alunos_met2 = 0; % Contador de sequencias feitas por alunos.
cont_matlab_met2 = 0; % Contador de sequencdias feitas pelo Matlab.

% Inic. vetor que armazena o resultado do teste de cada sequencia.
vect_result_met2 = zeros(1,tam_seq);

% Calculando o numero de digitos das sequencias (geralmente 100).
tam_seq = length(list_seq(1,:));

% Inicializando matriz que armazena as contagem das combinacoes possiveis.
mat_cont = zeros(tam_seq, 2^n_bits);

% Baseado em que todas combinacoes devem ter mesma quantidade de aparicoes.
cont_ideal = (100/n_bits)/2^n_bits;

margem = 0.5*cont_ideal; % Margem sobre o ideal para o teste de aceitacao.

% Loop que itera sob. todas as seq. realizando o teste de aleatoriedade.
for i = 1:num_seq
    % Varre a sequencia contagem o numero de aparicoes de cada combinacao.
    for j = 1:n_bits:tam_seq
        capture = bi2de(list_seq(i,j:j+n_bits-1), 'left-msb');
        mat_cont(i,capture+1) = mat_cont(i,capture+1)+1;
    end
%     figure(3);
%     bar(mat_cont(i,:),'k'); hold on;
%     plot([0 5], [cont_ideal+margem cont_ideal+margem], 'b');
%     plot([0 5], [cont_ideal-margem cont_ideal-margem], 'b'); 
%     plot([0 5], [cont_ideal cont_ideal], '--'); hold off;
%     pause
    % Analisa a contagem de cada combina��o e se estiver dentro da faixa
    % da quant. ideal + e - margem, incrementa a variavel do crit�rio de 
    % aleatoriedade.
    num_comb_dentro = 0;    
    for k = 1:length(mat_cont(i,:))
        if (mat_cont(i,k)<cont_ideal+margem)&&(mat_cont(i,k)>cont_ideal-margem)
            num_comb_dentro = num_comb_dentro + 1;
        end
    end
   
    % Conforme o n�mero de comb. dentro da faixa � visto a sequ�ncia como
    % gerada por aluno (n�o aleat�ria) ou por Matlab (aleat�ria).
    if num_comb_dentro <= 3
        cont_alunos_met2 = cont_alunos_met2 + 1;
        vect_result_met2(i) = 0;
    else
        cont_matlab_met2 = cont_matlab_met2 + 1;
        vect_result_met2(i) = 1;
    end     
end

figure(2);
stem(find(vect_result_met2==0), vect_result_met2(find(vect_result_met2==0)), 'r');
hold on;
stem(find(vect_result_met2>0), vect_result_met2(find(vect_result_met2>0)), 'g');
hold off;
legend(sprintf('Pseudo N�o-Aleat�rias (Total - %0.2d)',cont_alunos_met2), ...
       sprintf('Pseudo Aleat�rias (Total - %0.2d)', cont_matlab_met2));
str = [sprintf('Falso Positivo - %0.2d', length(find(vect_result_met2(1:35)==1)))...
       newline ...
       sprintf('Falso Negativo - %0.2d', length(find(vect_result_met2(36:70)==0)))];
annotation('textbox',[.14 .805 .1 .1],'String',str,'FitBoxToText','on');
%title('M�todo de contagem e teste de hip�tese baseado no Agrupamento de Bits');
axis([0 71 0 1.5]);

%==========================================================================

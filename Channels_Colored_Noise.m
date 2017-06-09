% =========================================================================
% Cap 9 - Teoria da Informação (pág. 277)
% 
% Canais com Ruído Colorido (CHANNELS WITH COLORED GAUSSIAN NOISE)
% 
% Objetivo: Encontrar Kx que maximiza a Capacidade do Canal.
% =========================================================================

clc;
clear all;
close all;

% =========================================================================

% Canal com pot. de ruído iguais e sem memória.
Kz_noNoise = [1 0 0; 0 1 0; 0 0 1]; 
% Canal com pot. de ruído diferentes e com memória.
Kz_MediumNoise = [1 0.3 0.2; 0.3 0.5 0.1; 0.2 0.1 2];
% Canal com pot. de ruído diferentes e com memória.
Kz_HighNoise = [1 0.3 0.2; 0.3 0.5 0.1; 0.2 0.1 0.9];

nP = 3; % Limite de potência no transmissor.

C = zeros(1,3);

for k = 1:3
    
    % A cada iteração para encontrar a capacidade troca-se o ruído.
    if k == 1 
        Kz = Kz_noNoise;
    else if k == 2
            Kz = Kz_MediumNoise;
         else
            Kz = Kz_HighNoise;
        end
    end
    
    [S,D] = eig(Kz); % Diagnolização da matriz, por uma matriz de similaridade.
    [Q,R] = qr(S); % Encontrar o congruente correspondente da matriz.
    
    lambda = [D(1,1) D(2,2) D(3,3)]; % Autovalores (valores da matriz diagnolizada).
    gamma = 0:0.01:10; % Faixa de gamma (nível de água) valores que será varrida.
    
    Aii = zeros(length(gamma),3); % Diagonal da matriz A: A = Q'*Kx*Q.
    sum_Aii = zeros(1,length(gamma)); % Vetor com o valores dos somatórios de Aii.
    
    % Parte do código que aplica o cálculo do máximo: Aii = (Gamma - lambda)^+
    for i = 1:length(gamma)
        temp = 0;
        for j = 1:3
            Aii(i,j) = (gamma(i) - lambda(j));
        end
        sum_Aii(i) = sum(Aii(i,:));
        if roundn(sum_Aii(i),-1) == nP
            index_max = i;
            break
        end
    end
    
    % Dúvida, sobre como encontrar o resto da matriz A, porque só
    % conhecemos a diagonal de A:
    A = [Aii(index_max,1) 0 0;0 Aii(index_max,2) 0;0 0 Aii(index_max,3)];
    
    % Transformada a matriz A em Kx através das matrizes de congruencia.
    Kx = Q*A*Q';
    
    % Variável para testar se o traço da matriz respeita o limite de
    % potência.
    testTrace = trace(Kx);
    
    % Cálculo da capacidade "conjunta" desse canal com Kx que fornece a
    % máxima I(X1,X2,X3;YI,Y2,Y3).
    C(k) = (1/2)*log2(det(Kx+Kz)/det(Kz)) 
end

% Plota o gráfico da capacidade "conjunta" desse canal para cada tipo de
% ruído.
figure(1);
plot(linspace(1,3,3), C, 'ob', linspace(1,3,3), C, 'r');

title('Channels with Colored Gaussian Noise');
xlabel('Noise Memory (Kz): 1 -> None, 2 -> Medium, 3 -> High');
ylabel('Capacity Maximum');

% =========================================================================
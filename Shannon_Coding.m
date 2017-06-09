% =========================================================================
% CAP 7 - Teoria da Informação
% 
% Código de Shannon - Sequência Conjuntamente Tipico
% 
% Objetivo: Construir um codificador de Shannon.
% =========================================================================

clc;
clear all;
close all;

% =========================================================================
% =========================================================================
%                           Transmissor

R = 0.5; % Taxa de bits (Taxa de informação)
alfa = 0.1; % Probabilidade de erro do canal BSC
n = 5; % Tamanho da sequência do bloco enviado

nComb = ceil(2^(n*R)) % Número de palavras códigos criadas

pX = [0.5 0.5];

C = randsrc(nComb,n,[0 1; pX(1) pX(2)]) % Palavras códigos geradas

w = randsrc(1,1,[1:1:nComb]);

x = C(w,:)

% =========================================================================
% =========================================================================
%                                Canal

y = bsc(x,alfa)

% =========================================================================
% =========================================================================
%                              Receptor                 
% 
% Conjuntamente Típico Decodificador

epson = 0.2;

Hx = 1; Hy = 1;

HCond_YX = zeros(1,length(alfa));
pCond_YX = [1-alfa alfa; alfa 1-alfa];

for i = 1:2 % Cálculo da entropia condicional H(Y|X).
    temp = 0;    
    for j = 1:2
        temp = temp +(-pCond_YX(j,i)*log2(pCond_YX(j,i)));
    end
    HCond_YX = HCond_YX + (temp*pX(i));
end

for i = 1:n % Cálculo da entropia condicional H(Y|X).
    temp = 0;    
    for j = 1:2
        temp = temp +(-pCond_YX(j,i)*log2(pCond_YX(j,i)));
    end
    HCond_YX = HCond_YX + (temp*pX(i));
end

Hxy = Hx + HCond_YX
Ixy = Hy - HCond_YX

lowBound = Hxy - epson
topBound = Hxy + epson

pxy = zeros(1,n);

for i = 1:nComb
    
    for j = 1:n
        for k = 1:2
            temp = pCond_YX(C(i:j))*pX(k);
        end
        pxy(j) = pxy(j)*pCond_YX(k,k)*pX(k);
    end
    In(i) = -(1/n)*log2(pxy(i));
end

k = 


% =========================================================================
% =========================================================================

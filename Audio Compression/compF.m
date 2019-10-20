function [som_comp,fs,CF,SNR] = compF(nome,nivel,tam_nivel,taxa,hab_graf)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% nivel de profundidade da WPD = 8 gerando 2^8 = 256 folhas
% tamanho das sub-bandas da WPD = 8 gerando 2^8 = 256 coeficientes
% taxas = [ ]32 [ ] 48 [ ] 64 [ ] 96 [ ] 128 [ ] 160 [ ] 192 [kbps]
% habilitação dos gráficos = [ ] 0-Desabilita [ ] 1-Habilita
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

dwtmode('per');
familia_wavelet = 'db40';

%%%%%%%%% Cálculo das matrizes para estimativa mínimos quadrados %%%%%%%%%%

eixo_x = -2^(tam_nivel-1):2^(tam_nivel-1)-1;
eixo_x = (eixo_x - mean(eixo_x))/std(eixo_x);
G_init = (eixo_x').^(0:15);

G1  = G_init(:,2:16);   % pol. ordem 15 com full termos e s/ zero
G2  = G_init(:,2:14);   % pol. ordem 13 com full termos e s/ zero
G3  = G_init(:,2:12);   % pol. ordem 11 com full termos e s/ zero
G4  = G_init(:,2:10);   % pol. ordem 09 com full termos e s/ zero
G5  = G_init(:,2:2:16); % pol. ordem 15 com só termos impares e s/ zero
G6  = G_init(:,2:2:14); % pol. ordem 13 com só termos impares e s/ zero
G7  = G_init(:,2:2:12); % pol. ordem 11 com só termos impares e s/ zero
G8  = G_init(:,2:2:10); % pol. ordem 09 com só termos impares e s/ zero
G9  = G_init(:,2:2:8);  % pol. ordem 07 com só termos impares e s/ zero
G10 = G_init(:,2:2:6);  % pol. ordem 05 com só termos impares e s/ zero
G11 = G_init(:,2:2:4);  % pol. ordem 03 com só termos impares e s/ zero

n_param_vector = [15 13 11 9 8 7 6 5 4 3 2];

%%%%%%%%%%%%%%%%%% Leitura do arquivo de entrada .wav %%%%%%%%%%%%%%%%%%%%%

% Lendo o arquivo de áudio
[som,fs] = audioread(strcat('.\audios_originais\',nome,'.wav'));

% Mixando caso o som seja estereo
if size(som,2) == 2
    som = (som(:,1)+som(:,2))/2;
elseif size(som,2)~=1
    error('Erro na dimensao do arquivo wav, so aceita mono ou estereo');
end

% Calculando a margem e tamanho das janelas e duração total do áudio
margem = floor(sqrt(2^nivel*2^tam_nivel));
duracao = 2^nivel*2^tam_nivel-2*margem;
duracao_total = length(som);

% Deve-se garantir que "excesso" não seja nulo no último chunk
if mod(duracao_total,duracao) == 0
    som(end)=[];
end

% Normalizando o sinal de áudio de -1 a 1
fator_norm = max(abs(som));
som_norm = som/fator_norm;

%%%%%%%%%%%% Obtendo os valores do modelo de limiar de audição %%%%%%%%%%%%

limiar_modelo = threshold_model(nivel,fs);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% Construção do cabeçalho no arquivo de saída .cww %%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Criando o arquivo de áudio comprimido
fid = fopen(strcat('.\audios_comprimidos\',nome,'.cww'),'w');

% Escrevendo os dados de cabeçalho
fwrite(fid,[nivel tam_nivel taxa],'uint8');
fwrite(fid,fator_norm,'float32');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Criando variável que armanezara o valor do som comprimido
som_comp = [];

% Criando apontador para dizer até onde no som processamos
apontador = 1;

% Obtendo o áudio em chunks e comprimindo
while apontador <= duracao_total
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%% Segmentação do Áudio %%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    som_chunk = zeros(max(margem-apontador+1,0),1);
    som_chunk = [som_chunk;som_norm(max(apontador-margem,1):min(apontador-1+duracao+margem,duracao_total))];
    
    % Calcula o excesso que será preenchido para formar o ultimo chunk
    excesso = max(apontador-1+duracao-duracao_total,0);
    
    % Grava o valor de excesso no cabeçalho do arquivo de áudio
    fwrite(fid,excesso,'uint32');
    
    %%%%%%%%%%%%%%% Calculando o tamanho do buffer por frame %%%%%%%%%%%%%%%%%%
    buffer = floor((duracao-excesso)/fs*taxa*1000);
    
    % Atribuindo os bits consumidos "fixos" dos cabeçalhos dos frames
    % 32-excesso 10-limiar 4-G_index  tam_nivel-fator_peso
    buffer = buffer-(32+10+4+tam_nivel);  

    % Se necessário completa o chunk com zeros
    som_chunk = [som_chunk;zeros(duracao+2*margem-length(som_chunk),1)];
    apontador = min(apontador+duracao,duracao_total+1);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%% Decomposição em wavelets packet %%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    tree_wpd = wpdec(som_chunk, nivel, familia_wavelet);
    folhas = leaves(tree_wpd);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%% Processo de procura pelos parâmetros adequados %%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    estado = 0;           % estado inicial                 
    G_index = 7;          % ordem inicial
    ajuste_limiar = 0.1;  % limiar inicial
    fator_peso = 1;       % fator de peso inicial
    
    gasto_ind = tam_nivel;
    
    if mod(gasto_ind,8) ~= 0
        gasto_ind = gasto_ind+(8-mod(gasto_ind,8));
    end
    
    flag_fim = false;     
    
    % Obtendo o máximo de cada folha e deixando armazenado
    fls_max_vector = zeros(1,length(folhas));
    for k = 1:length(folhas)
        i = folhas(k);
        y = wpcoef(tree_wpd,i);
        fls_max_vector(k) = max(abs(y));
    end            

    while flag_fim == false
                
        if ismember(estado,[-2 0 +2 +3 +4])    
            
            if apontador == 1
                % Cabeçalho (8-(nivel,tam_nivel e taxa) 32-fator_norm)
                buffer_A = buffer-(3*8+32); 
            else
                buffer_A = buffer;
            end      
            
            fora = [];
            dentro = [];  
            
            limiar = limiar_modelo*ajuste_limiar; 
            
            for k = 1:length(folhas)
                if fls_max_vector(k) < limiar(k)
                    fora = [fora k-1];
                else
                    dentro = [dentro k-1];
                end
            end            
            
            buffer_A = buffer_A-(nivel+1+1);
            
            if length(fora) >= length(dentro)
                if not(isempty(dentro))
                    dentro_envio = byte_cast(dentro,nivel,8);
                    buffer_A = buffer_A-(length(dentro_envio)*8);
                end
            else
                if not(isempty(fora))
                    fora_envio = byte_cast(fora,nivel,8);
                    buffer_A = buffer_A-(length(fora_envio)*8);
                end
            end
            
            n_coef = zeros(1,length(folhas));
            
            for k = dentro
                i = folhas(k+1);
                y = wpcoef(tree_wpd,i);
                pos = find(y<=-limiar(k+1) | y>=limiar(k+1));
                n_coef(k+1) = length(pos);
            end
        end
        
        buffer_B = buffer_A;
        
        for k = dentro    
            n_param = n_param_vector(G_index); 
            
            gasto_posicoes = n_coef(k+1)*tam_nivel;             
            if mod(gasto_posicoes,8) ~= 0
                gasto_posicoes = gasto_posicoes+(8-mod(gasto_posicoes,8));
            end                        
            
            if (n_coef(k+1)-fator_peso)*16 <= ((n_param)*32+gasto_ind+2)
                buffer_B = buffer_B-(n_coef(k+1)*16+gasto_posicoes+gasto_ind+8);
            else
                buffer_B = buffer_B-(n_param*32+gasto_posicoes+2*gasto_ind+8);
            end
        end        

        [estado,G_index,ajuste_limiar,fator_peso,flag_fim] = ...
        state_machine(estado,buffer_B,flag_fim,G_index,ajuste_limiar,fator_peso,tam_nivel);  
    
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%% Codificando o Chunk utilizando %%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%% os parâmetros encontrados %%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    disp(num2str([estado,n_param_vector(G_index),ajuste_limiar,fator_peso,flag_fim],...
    'Estado = %.0f  #Polinômios = %.0f  Limiar = %.3f  Fator Peso = %.0f  flag fim = %.0f'));
    
    % Escrevendo no cabeçalho do frame os parâmetros finais encontrados:
    fwrite(fid,[de2bi(uint16(ajuste_limiar*1000),10) de2bi(G_index,4) de2bi(fator_peso-1,tam_nivel)],'ubit1');
    
    switch G_index
        case 1;  G = G1;  case 2;  G = G2; case 3;  G = G3;
        case 4;  G = G4;  case 5;  G = G5; case 6;  G = G6;
        case 7;  G = G7;  case 8;  G = G8; case 9;  G = G9;
        case 10; G = G10; case 11; G = G11;
        otherwise; error('erro na programacao dos estados');
    end
    
    limiar = limiar_modelo*ajuste_limiar;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%% Eliminação de sub-bandas %%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % Criando árvore WPD com valores nulos p/ transformada inversa.
    tree_rec = wpdec(zeros(duracao+2*margem,1),nivel,familia_wavelet);
    
    fora = [];
    dentro = [];
    
    % Elimando as folhas c/ máximo valor abaixo do limiar
    for k = 1:length(folhas)
        i = folhas(k);
        y = wpcoef(tree_wpd,i);
        max_f = max(abs(y));
        if max_f < limiar(k)
            fora = [fora k-1];
        else
            dentro = [dentro k-1];
        end
    end
    
    fwrite(fid,de2bi(length(dentro),nivel+1),'ubit1');
    
    if length(fora) >= length(dentro)
        % envia ao cabeçalho somente as folhas que ficam dentro
        if not(isempty(dentro))
            fwrite(fid,byte_cast(dentro,nivel,8),'uint8');
        end
    else
        % envia ao cabeçalho somente as folhas que ficam fora
        if not(isempty(fora))
            fwrite(fid,byte_cast(fora,nivel,8),'uint8');
        end
    end
    
    % Varrendo as folhas sobreviventes para saber qual codificar ou não
    for k = dentro      
        
        i = folhas(k+1);
        y = wpcoef(tree_wpd,i);
        pos = find(y<=-limiar(k+1) | y>=limiar(k+1));
        n_coef = length(pos);
        n_param = size(G,2);
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%% Codificando ou não codificando as sub-bandas %%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        % Dependendo do n° de coef. codifica ou não as sub-bandas:
        if (n_coef-fator_peso)*16 <= ((n_param)*32+gasto_ind+2)
            
            flag_cod = 0;
            % Atribuindo valores dos coeficientes
            coefs = y(pos);
            % Encontrando valor máximo dos coeficientes para quantização
            max_coef = max(abs(coefs));
            if max_coef >= 30
                error('Valor máximo (fator de escala) não esperado');
            end
            % Valor mínimo para quantização é o limiar da sub-banda
            min_quant = limiar(k+1);
            
            % Quantizando o valor máximo (fator de escala) em 7 bits.
            passo_quant = (30-min_quant)/2^7;
            max_quant = uint8(ceil((max_coef-min_quant)/passo_quant));
            
            if max_quant >= 128
                max_quant = 127;
            end
            
            % Recuperando valor máximo quantizado
            max_aprox = double(max_quant)*passo_quant + min_quant;
            
            % Quantizando os coeficientes com 16 bits
            coefs_quant = zeros(n_coef,1);
            passo_quant_coefs = (max_aprox-min_quant)/2^(16-1);
            
            for p = 1:n_coef
                if coefs(p) >= 0
                    coefs_quant(p) = int16((coefs(p)-min_quant)/passo_quant_coefs);
                else
                    coefs_quant(p) = int16((coefs(p)+min_quant)/passo_quant_coefs);
                end
            end
           
            fwrite(fid,[flag_cod de2bi(max_quant,7)],'ubit1');
            fwrite(fid,coefs_quant,'uint16');
            fwrite(fid,byte_cast(n_coef-1,tam_nivel,8),'uint8');
            fwrite(fid,byte_cast(pos-1,tam_nivel,8),'uint8');
            
        else
            
            flag_cod = 1;
            
            [z,ind] = sort(y,'ascend');
            
            N1_estouro_flag = 0;
            N2_estouro_flag = 0;
            
            N1 = length(find(z<=-limiar(k+1)));
            N2 = length(find(z>=limiar(k+1)));

            tam_subbanda = length(z);
            
            % Matriz de pesos
            pesos = ones(tam_subbanda,1);
            
%             if n_param >= 9
%                 if N1>0
%                     pesos(N1)=1;
%                     pesos(1)=10;    
%                 end
%                 if N2>0
%                     pesos(end-N2+1)=1;
%                     pesos(end)=10;
%                 end
%             end
            
            P = diag(pesos);
            
            % Cálculo das estimativas mínimos quadrados
            theta = single((P*G)\(P*z));
            w = (G*theta);

            if N1 > 0
                ind_abaixo = 1:N1;
            else
                ind_abaixo = [];
            end
            if N2 > 0
                ind_acima = (tam_subbanda-N2+1):tam_subbanda;
            else
                ind_acima = [];
            end    
            
            if N1 == 2^tam_nivel
                N1_estouro_flag = 1;
            end            
            if N2 == 2^tam_nivel
                N2_estouro_flag = 1;
            end
            
            posicoes = [ind(ind_abaixo)-1; ind(ind_acima)-1];
            
            fwrite(fid,[flag_cod N1_estouro_flag N2_estouro_flag],'uint8');
            fwrite(fid,theta,'float32');
            fwrite(fid,byte_cast([N1 N2],tam_nivel,8),'uint8');
            fwrite(fid,byte_cast(posicoes,tam_nivel,8),'uint8');
            
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%% Reconstrução do sinal de saída %%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        if flag_cod == 0
            coef_desquant = zeros(n_coef,1);
            for p = 1:n_coef
                if coefs_quant(p) >= 0
                    coef_desquant(p) = double(coefs_quant(p))*passo_quant_coefs + min_quant;
                else
                    coef_desquant(p) = double(coefs_quant(p))*passo_quant_coefs - min_quant;
                end
            end
            r = zeros(2^tam_nivel,1);
            r(pos,1) = coef_desquant;
        else
            ind_nulos = (N1+1):(tam_subbanda-N2-1);
            w(ind_nulos) = 0;          
            for j = 1:length(ind)
                r(ind(j),1) = w(j);
            end
        end
        tree_rec = write(tree_rec,'cfs',i,r);
    end
    
    s1 = wprec(tree_rec);
    s2 = s1(margem+1:end-margem-excesso);
    som_comp = [som_comp; s2*fator_norm];
    
%     fwrite(fid,zeros(buffer,1),'ubit1');
end

erro = (som-som_comp);
fclose(fid);

%%%%%%%%%%%%%%%%%%% Parâmetros de avaliação de desempenho %%%%%%%%%%%%%%%%%

SNR = snr(som,erro);

D1 = dir(strcat('.\audios_comprimidos\',nome,'.cww'));
D2 = dir(strcat('.\audios_originais\',nome,'.wav'));
CF = D2.bytes/D1.bytes;

% PSNR_dB = 10*log10((max(x)^2)/erro_mse);

%%%%%%%%%%%%%%%%%%%%%%% Gráficos de saída %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if hab_graf
    sound(som_comp,fs);
    tempo = (1:duracao_total)/fs;
    
    subplot(2,1,1);
    plot(tempo,som); hold on; grid on;
    plot(tempo,som_comp,'r'); hold off;
    axis([0 tempo(end) -max(abs(som)) max(abs(som))]);
    
    subplot(2,1,2);
    plot(tempo,erro,'k'); grid on;
    axis([0 tempo(end) -max(abs(som)) max(abs(som))]);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
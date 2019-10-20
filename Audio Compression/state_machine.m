function [estado,index_G,ajuste_limiar,fator_peso,flag_fim] = ...
  state_machine(estado,buffer,flag_fim,index_G,ajuste_limiar,fator_peso,tam_nivel)
%==========================================================================
% Função recebe o estado atual e toma a decisão da escolha do próximo 
% estado através da matriz de transição: 
%==========================================================================

ajuste_limiar = round(ajuste_limiar,3);

if estado == 0
    if buffer > 0
        estado = +1;
        index_G = index_G - 1;
    elseif buffer <= 0
        estado = -1;
        index_G = index_G + 1;
    end
    
%==========================================================================
% Ajustando a ordem polinomial com buffer aproximado (qualidade ou compressão)

elseif estado == +1 % (qualidade)
    if buffer > 0
        if index_G > 1
            index_G = index_G - 1;
        else
            estado = +3;
            ajuste_limiar = ajuste_limiar - 0.01;
        end
    else
        index_G = index_G + 1;
        estado = +3;
        ajuste_limiar = ajuste_limiar - 0.01;
    end
    
elseif estado == -1 % (compressão)
    if buffer < 0
        if index_G < 11
            index_G = index_G + 1;
        else
            estado = -2;
            ajuste_limiar = ajuste_limiar + 0.1;
        end
    else
        estado = +3;
        ajuste_limiar = ajuste_limiar - 0.01;
    end     
    
%==========================================================================
% aumenta o "valor inicial" do limiar em incrementos de 0.1 (compressão)

elseif estado == -2
    if buffer < 0
        ajuste_limiar = ajuste_limiar + 0.1;
    else
        estado = +3;
        ajuste_limiar = ajuste_limiar - 0.01;
    end
    
%==========================================================================
% diminui o valor inicial em incrementos de -0.01 (qualidade)

elseif estado == +3 
    if buffer > 0
        if ajuste_limiar > 0.01   
            ajuste_limiar = ajuste_limiar - 0.01;
        else
            estado = +5;
            fator_peso = fator_peso + 1;
        end
    else
        ajuste_limiar = ajuste_limiar + 0.01;
        estado = +4;
        ajuste_limiar = ajuste_limiar - 0.001;
    end   
    
%==========================================================================
% diminui o valor inicial em incrementos de -0.001 (qualidade)    

elseif estado == +4 
    if buffer > 0
        ajuste_limiar = ajuste_limiar - 0.001;
    else
        ajuste_limiar = ajuste_limiar + 0.001;
        flag_fim = true;
    end

%==========================================================================
% aumenta o fator_peso em incrementos de -17 (qualidade)    

elseif estado == +5 
    if buffer > 0
        if fator_peso < (2^tam_nivel)
            fator_peso = fator_peso + 1;
        else
            estado = +4;
            ajuste_limiar = ajuste_limiar - 0.001;
        end
    else
        fator_peso = fator_peso - 1;
        estado = +4;
        ajuste_limiar = ajuste_limiar - 0.001;
    end
end

return
%==========================================================================

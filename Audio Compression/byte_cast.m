function seq_saida = byte_cast(seq_entrada,nbits_entrada,nbits_saida)

if max(seq_entrada)>2^nbits_entrada-1 || min(seq_entrada)<0
    disp('erro')
end

j = 1;
apont = 1;
seq_saida(1) = 0;

for i=1:length(seq_entrada)    
    resto = nbits_entrada;
    while resto > 0
        bits_envio = min((nbits_saida-apont+1),resto);
        resto = max(resto-(nbits_saida-apont+1),0);
        seq_saida(j) = seq_saida(j)+bitshift(bitand(seq_entrada(i),...
                       2^(bits_envio+resto)-1),...
                       nbits_saida-apont+1-bits_envio-resto);
        apont = mod(apont-1+bits_envio,nbits_saida)+1;
        if apont == 1 && (i < length(seq_entrada) || resto > 0)
            j = j+1;
            seq_saida(j) = 0;
        end
    end
end

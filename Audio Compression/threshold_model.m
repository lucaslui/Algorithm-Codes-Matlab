function limiar = threshold_model(nivel,fs)
%==========================================================================
% Função recebe o nível de profundidade da transformada wavelet e a 
% frequencia de amostragem do sinal de áudio e retorna o modelo do limiar
% de audição padronizado pela ISO226 ajustado as sub-bandas wavelets.
%==========================================================================

num_fls = 2^nivel;

% valores da ISO226 do limiar auditivo
[spl_iso, freq_iso] = iso226(10);

log_f = log2([freq_iso(19:end) 14200 22050]);
log_f = [4 10 log_f];
db_amp = [20 0 spl_iso(19:end) 9 20];

a = zeros(length(db_amp)-1,1);
b = zeros(length(db_amp)-1,1);

for i = 1:length(log_f)-1
    if log_f(i) > log_f(i+1)
        erro('erro curva do limiar')
    end
    a(i) = (db_amp(i+1)-db_amp(i))/(log_f(i+1)-log_f(i));
    b(i) = db_amp(i)-log_f(i)*a(i);
end

limiar = zeros(num_fls,1);

for i = 1:num_fls
    freq = (i-0.5)*fs/(2*num_fls);
    log_freq = log2(freq);
    j = length(a);
    while (log_freq<log_f(j)) && (j>1)
        j = j-1;
    end
    db_amplitude = a(j)*log_freq+b(j);
    limiar(i) = 10^(db_amplitude/20);
end
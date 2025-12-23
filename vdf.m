function ni = vdf(n, T, AD)
    % n - общая концентрация (скаляр)
    % T - температура (скаляр)
    
    beta = 1 ./ (AD.k * T);
    expE = exp(-AD.e1234 * beta);
    
    Zv = sum(AD.stw .* expE); 
    ni = (n / Zv) .* (AD.stw .* expE);
end
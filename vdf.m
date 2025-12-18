function ni = vdf(n, T, AD)
% функция расчёта колебательных распределений молекул СH4 при температуре T [К]

expE = exp(-AD.e1234 * (AD.k * T).^(-1));
Zv = sum(AD.stw .* expE); 
ni = (n ./ Zv) .* (AD.stw .* expE);

end

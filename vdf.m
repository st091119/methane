function ni = vdf(n, T, AD)
% функция расчёта колебательных распределений молекул СH4
% входные параметры содержат вектора-столбцы переменных

[ii,jj,kk, ll] = ind2sub(AD.lch4, AD.ind_1d_in_4d);
i_state = ii-1; %(symmetrical)
j_state = jj-1; %(twisting)
k_state = kk-1; %(antisymmetrical)
l_state = ll-1; %(scissoring)

% распределение Больцмана
E2 = i_state*AD.e1000 * (AD.k*T').^(-1);
E3 = j_state*AD.e0100 * (AD.k*T').^(-1);
E4 = k_state*AD.e0010 * (AD.k*T').^(-1);
E5 = l_state*AD.e0001 * (AD.k*T').^(-1);

expE = exp(-E2-E3-E4-E5);
Zv = sum(AD.stw .* expE); 
ni = (n'./Zv).*(AD.stw.*expE);

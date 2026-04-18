function ni = vdf_multitemperature(n, T1, T2, T3, T4, AD)
%VDF_MULTITEMPERATURE  Заселённости CH4 по мультитемпературной формуле (аналог (19)).
%
%   ni = vdf_multitemperature(n, T1, T2, T3, T4, AD)
%
%   n_{i1,i2,i3,i4} ∝ s_{i1,i2,i3,i4} * exp( -Σ_m i_m * ε_m / (k_B T_m) ),
%   где i_m — квантовые числа из AD.inds (0-based), ε_m — энергия одного кванта
%   моды m от основного состояния: ε_m = E(…,1_m,…) - E(0000) из input_data.
%
%   Вход:
%     n      — масштаб плотности [м^-3] (как в vdf.m)
%     T1..T4 — колебательные температуры мод ν1…ν4 [К]
%     AD     — структура из input_data (stw, inds, k, e0000, e1000, …)
%
%   Выход:
%     ni — плотность на уровне [м^-3], сумма по уровням = n.

Tm = [T1, T2, T3, T4];
Tm = max(Tm(:).', 1e-300);

kB = AD.k;
e0 = AD.e0000;
eps = [AD.e1000 - e0, AD.e0100 - e0, AD.e0010 - e0, AD.e0001 - e0];

I0 = AD.inds(:, 1);
J0 = AD.inds(:, 2);
K0 = AD.inds(:, 3);
L0 = AD.inds(:, 4);

xi = -(I0 * eps(1) / (kB * Tm(1)) + J0 * eps(2) / (kB * Tm(2)) + ...
    K0 * eps(3) / (kB * Tm(3)) + L0 * eps(4) / (kB * Tm(4)));

expE = exp(xi);
Zv = sum(AD.stw .* expE);
ni = (n ./ Zv) .* (AD.stw .* expE);

end

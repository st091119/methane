function ni = vdf_mt(n, T, Tv2, Tv3, Tv4, AD)
%VDF_MT  Распределение Больцмана CH4 при T для моды 1 и Tv2,Tv3,Tv4 для мод 2–4.
%   ni = vdf_mt(n, T, Tv2, Tv3, Tv4, AD)
%   Частицы на уровень пропорциональны AD.stw .* exp(-E1/(kT) - E2/(kTv2) - ...).

hc = AD.h * AD.c;
I0 = AD.inds(:, 1);
J0 = AD.inds(:, 2);
K0 = AD.inds(:, 3);
L0 = AD.inds(:, 4);
d = AD.d;

E1 = hc .* AD.omega(1) .* (I0 + d(1)/2);
E2 = hc .* AD.omega(2) .* (J0 + d(2)/2);
E3 = hc .* AD.omega(3) .* (K0 + d(3)/2);
E4 = hc .* AD.omega(4) .* (L0 + d(4)/2);

xi = -E1 ./ (AD.k * T) - E2 ./ (AD.k * Tv2) - E3 ./ (AD.k * Tv3) - E4 ./ (AD.k * Tv4);
expE = exp(xi);
Zv = sum(AD.stw .* expE);
ni = (n ./ Zv) .* (AD.stw .* expE);

end

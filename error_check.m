function d1 = error_check(Y0, Y1, AD)
% проверка выполнения закона сохранения энергии: E = const

% параметры газа в начальный момент времени и после
Y0 = Y0(:);
if numel(Y0) < 2
	error('error_check: InvalidY0','Y0 must contain [T0 Tv0] or [T0 Tv2 Tv3 Tv4].');
end
T0 = Y0(1);

Y1 = Y1(:);
if numel(Y1) < 2
	error('error_check: InvalidY1','Y1 must contain [T Tv] or [T Tv2 Tv3 Tv4].');
end
T = Y1(1);

use_mt = (numel(Y0) >= 4 && numel(Y1) >= 4);
if use_mt
    Tv2_0 = Y0(2); Tv3_0 = Y0(3); Tv4_0 = Y0(4);
    Tv2 = Y1(2);   Tv3 = Y1(3);   Tv4 = Y1(4);
else
    Tv_0 = Y0(2);
    Tv = Y1(2);
end

% безразмерные колебательные распределения молекул
n0 = 1 %AD.n0;
n = 1 %AD.p0 / (AD.k * T); % плотность частиц [м^-3]

if use_mt
    ni = vdf_mt(n, T, Tv2, Tv3, Tv4, AD);
    ni_0 = vdf_mt(n0, T0, Tv2_0, Tv3_0, Tv4_0, AD);
else
    ni = vdf(n, Tv, AD);
    ni_0 = vdf(n0, Tv_0, AD);
end

% колебательная энергия
E = AD.e1234; % [Дж]
e_v = sum(ni .* E);
e_v0 = sum(ni_0 .* E);

% поступательная энергия смеси
e_tr = 1.5 * AD.k * T * n;
e_tr0 = 1.5 * AD.k * T0 * n0;

% вращательная энергия CH4 (модель жесткого ротора)
e_rot = 1.5 * AD.k * T * n;
e_rot0 = 1.5 * AD.k * T0 * n0;

% внутренняя энергия газа [Дж/м^3]
E = e_tr + e_rot + e_v;
E0 = e_tr0 + e_rot0 + e_v0;

% невязка
u1 = E0 - E;

% относительная невязка
d1 = max(abs(u1) / E0);
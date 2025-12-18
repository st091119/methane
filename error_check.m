function d1 = error_check(Y0, Y1, AD)
% проверка выполнения закона сохранения энергии: E = const

% параметры газа в начальный момент времени и после
Y0 = Y0(:);
if numel(Y0) < 2
	error('error_check: InvalidY0','Y0 must contain [T0 Tv0].');
end
[T0, Tv_0] = deal(Y0(1), Y0(2));

Y1 = Y1(:);
if numel(Y1) < 2
	error('error_check: InvalidY1','Y1 must contain [T Tv].');
end
[T, Tv] = deal(Y1(1), Y1(2));

% безразмерные колебательные распределения молекул
n0 = 1 %AD.n0;
n = 1 %AD.p0 / (AD.k * T); % плотность частиц [м^-3]

ni = vdf(n, Tv, AD);
ni_0 = vdf(n0, Tv_0, AD);

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
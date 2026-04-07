function kup = kvt2_ssh(T, AD)
% KVT2_SSH коэффициенты скорости VT-обмена v2 -> v2+1 (мода nu_2) при CH4-CH4.
% Теория SSH по постановке как для CO2-CO2 (Adamovich / da Silva), с константами CH4.
%
%   T  — температура [К]
%   AD — структура из input_data (+ поля, добавляемые в main при необходимости)
%
% Если заданы AD.indjp (длина N, NaN где нет соседа), используются они; иначе соседи
% строятся по AD.inds: квантовое число второй моды J -> J+1 при неизменных I,K,L.
%
% Выход: kup [м^3/с] — коэффициенты скорости для переходов с заданным партнёром.

N = numel(AD.e1234);

if isfield(AD, 'indjp') && numel(AD.indjp) == N)
    indjp = AD.indjp(:);
else
    indjp = build_indjp_vt2(AD);
end

ip = find(~isnan(indjp));
ip1 = indjp(ip);

% квантовое число моды 2 (0-based), как «st_before» в исходной CO2-формуле
st_before = AD.inds(ip, 2);

% газодинамический радиус взаимодействия [м]
R0 = AD.r0;

% обратный радиус взаимодействия [м^-1]
alpha = 17.5 ./ R0;

% приведённая масса сталкивающихся молекул [кг] (идентичные CH4)
mu = 0.5 * AD.m;

delta_E = AD.e1234(ip) - AD.e1234(ip1);

% потенциал взаимодействия: одноквантовые переходы
A = 0.5;

% приведённая масса осциллятора по моде 2 [кг] — как в probabilties_fho_vt.m
M = 1.6734911e-27;

% частота колебаний моды nu_2 [с^-1]: omega в см^-1, c в см/с (input_data)
om_e_cm = AD.omega(2);
nu = om_e_cm * AD.c;

al = 4 * pi^2 * M * nu / AD.h;
V = -A * sqrt(st_before + 0.5 - 0.5 * sign(delta_E)) * alpha ./ ...
    sqrt(2 * al);

p = zeros(N, 1);

% нерезонансный обмен
sigma = 3 * (2 * pi^4 * delta_E.^2 .* (mu ./ (alpha.^2 * AD.h^2 * AD.k * T))).^(1/3) + ...
    delta_E ./ (2 * AD.k * T);

p(ip) = V.^2 * 0.394 .* (8 * pi^3 * delta_E .* (mu ./ (alpha * AD.h).^2)).^2 .* ...
    sigma.^1.5 .* exp(-sigma) ./ (1 - exp(-2/3 * sigma));

checker = p > 1;
p(checker) = 1;

% частота столкновений [м^3/с]
Z = R0.^2 ./ sqrt(mu) * sqrt(8 * pi * AD.k * T);
kup = Z .* p;

end

% --- локальная функция --------------------------------------------------------

function indjp = build_indjp_vt2(AD)
% Индекс партнёра при VT2+: (I,J,K,L) -> (I,J+1,K,L), только если уровень есть в AD.inds.

N = size(AD.inds, 1);
indjp = nan(N, 1);
Jmx = AD.lch4(2) - 1; % максимальное 0-based J

for r = 1:N
    v = AD.inds(r, :);
    if v(2) >= Jmx
        continue;
    end
    tgt = v + [0, 1, 0, 0];
    m = find(all(AD.inds == tgt, 2), 1);
    if ~isempty(m)
        indjp(r) = m;
    end
end

end

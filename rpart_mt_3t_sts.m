function dy = rpart_mt_3t_sts(~, y, AD)
%Трехтемпературная CH4-модель: y = [T_b; T13_b; T24_b].


if ~exist('rates_fho_vt', 'file')
    addpath(fullfile(fileparts(mfilename('fullpath')), 'fho_model'));
end

AD = ch4_relax_topology(AD);

T_b = max(y(1), 1e-9);
T13_b = max(y(2), 1e-9);
T24_b = max(y(3), 1e-9);

temp = T_b * AD.T0;
t13 = T13_b * AD.T0;
t24 = T24_b * AD.T0;

kT0 = AD.k * AD.T0;
N = numel(AD.e1234);

I0 = AD.inds(:, 1);
J0 = AD.inds(:, 2);
K0 = AD.inds(:, 3);
L0 = AD.inds(:, 4);

eps1 = AD.e1000 - AD.e0000;
eps2 = AD.e0100 - AD.e0000;
eps3 = AD.e0010 - AD.e0000;
eps4 = AD.e0001 - AD.e0000;

E13 = I0 * eps1 + K0 * eps3;
E24 = J0 * eps2 + L0 * eps4;

% Мультитемпературное распределение: T13 для (1,3), T24 для (2,4).
xi = -(E13 ./ (AD.k * t13) + E24 ./ (AD.k * t24));
f = AD.stw(:) .* exp(xi);
Zv = sum(f);
nco2i_b = f ./ Zv;

if ~isfield(AD, 'fho_steric_vv_34s')
    AD.fho_steric_vv_34s = 1.0;
end
if ~isfield(AD, 'fho_steric_vv_34d')
    AD.fho_steric_vv_34d = 0.06;
end
if ~isfield(AD, 'fho_gamma_vv34d')
    AD.fho_gamma_vv34d = 1.0;
end
if ~isfield(AD, 'fho_steric_vv_34_4d')
    AD.fho_steric_vv_34_4d = 1.0;
end

n_scale = AD.n0 * AD.tau;
sk0 = [0, 0, 0, 0];

Jmx = AD.lch4(2) - 1;
Lmx = AD.lch4(4) - 1;

k2_tab = zeros(Jmx + 1, 1);
for Jq = 1:Jmx
    k2_tab(Jq) = rates_fho_vt(temp, Jq, Jq - 1, 2, AD.fho_steric2, AD.fho_alpha, AD.fho_e_m, 'n_steps', 8000);
end

k4_tab = zeros(Lmx + 1, 1);
for Lq = 1:Lmx
    k4_tab(Lq) = rates_fho_vt(temp, Lq, Lq - 1, 4, AD.fho_steric4, AD.fho_alpha, AD.fho_e_m, 'n_steps', 8000);
end

RVT2 = zeros(N, 1);
src = find(~isnan(AD.indj_vt2));
for ii = 1:numel(src)
    r = src(ii);
    dst = AD.indj_vt2(r);
    Jq = AD.inds(r, 2);
    kf = k2_tab(Jq + 1) * n_scale;      % forward: dst(J+1) -> r(J)
    kr = k_vt_bwd(kf, dst, r, temp, AD); % backward: r(J) -> dst(J+1)
    RVT2(r) = RVT2(r) + nco2i_b(dst) * kf - nco2i_b(r) * kr;
    RVT2(dst) = RVT2(dst) + nco2i_b(r) * kr - nco2i_b(dst) * kf;
end

RVT4 = zeros(N, 1);
src = find(~isnan(AD.indl_vt4));
for ii = 1:numel(src)
    r = src(ii);
    dst = AD.indl_vt4(r);
    Lq = AD.inds(r, 4);
    kf = k4_tab(Lq + 1) * n_scale;      % forward: dst(L+1) -> r(L)
    kr = k_vt_bwd(kf, dst, r, temp, AD); % backward: r(L) -> dst(L+1)
    RVT4(r) = RVT4(r) + nco2i_b(dst) * kf - nco2i_b(r) * kr;
    RVT4(dst) = RVT4(dst) + nco2i_b(r) * kr - nco2i_b(dst) * kf;
end

RVV = zeros(N, 1);
src = find(~isnan(AD.indvv34));
kvv_cache = containers.Map('KeyType', 'char', 'ValueType', 'double');
for ii = 1:numel(src)
    r = src(ii);
    dst = AD.indvv34(r);
    si = AD.inds(r, :);
    sf = AD.inds(dst, :);
    key = sprintf('%d_%d', si(3), si(4));
    if isKey(kvv_cache, key)
        kv = kvv_cache(key);
    else
        kv = rates_fho_vv(temp, si, sf, sk0, sk0, AD.fho_steric_vv_34s, AD.fho_alpha, AD.fho_e_m, 'n_steps', 6000);
        kvv_cache(key) = kv;
    end
    kf = kv * n_scale;
    kr = k_vt_bwd(kf, r, dst, temp, AD);
    RVV(r) = RVV(r) + nco2i_b(dst) * kr - nco2i_b(r) * kf;
    RVV(dst) = RVV(dst) + nco2i_b(r) * kf - nco2i_b(dst) * kr;
end

RVV34d = zeros(N, 1);
src = find(~isnan(AD.indvv34d));
kvv34d_cache = containers.Map('KeyType', 'char', 'ValueType', 'double');
for ii = 1:numel(src)
    r = src(ii);
    dst = AD.indvv34d(r);
    si = AD.inds(r, :);
    sf = AD.inds(dst, :);
    key = sprintf('%d_%d', si(3), si(4));
    if isKey(kvv34d_cache, key)
        kv = kvv34d_cache(key);
    else
        kv = rates_fho_vv(temp, si, sf, sk0, sk0, AD.fho_steric_vv_34d, AD.fho_alpha, AD.fho_e_m, ...
            'n_steps', 6000, 'gamma', AD.fho_gamma_vv34d);
        kvv34d_cache(key) = kv;
    end
    kf = kv * n_scale;
    kr = k_vt_bwd(kf, r, dst, temp, AD);
    RVV34d(r) = RVV34d(r) + nco2i_b(dst) * kr - nco2i_b(r) * kf;
    RVV34d(dst) = RVV34d(dst) + nco2i_b(r) * kf - nco2i_b(dst) * kr;
end

[RVV34_4d, RVV34_4d_partner_E24] = ch4_vv34_4d_source(temp, nco2i_b, AD, n_scale, 6000);

% Веса инвариантов (как в формуле: i2*e0100+i4*e0001 и аналогично для 1,3).
W13 = E13 / kT0;
W24 = E24 / kT0;

% Явное разложение релаксационных членов по каналам.
R13_vv34 = sum(W13 .* RVV);
R13_vv34d = sum(W13 .* RVV34d);
R13_vv34_4d = sum(W13 .* RVV34_4d);
R13 = R13_vv34 + R13_vv34d + R13_vv34_4d;

R24_vt2 = sum(W24 .* RVT2);
R24_vt4 = sum(W24 .* RVT4);
R24_vv34 = sum(W24 .* RVV);
R24_vv34d = sum(W24 .* RVV34d);
R24_vv34_4d = sum(W24 .* RVV34_4d) + RVV34_4d_partner_E24 / kT0;
R24 = R24_vt2 + R24_vt4 + R24_vv34 + R24_vv34d + R24_vv34_4d;

% Аналитические производные средних энергий (без конечных разностей).
E13m = sum(nco2i_b .* E13);
E24m = sum(nco2i_b .* E24);
E13sq_m = sum(nco2i_b .* (E13.^2));
E24sq_m = sum(nco2i_b .* (E24.^2));
E13E24_m = sum(nco2i_b .* (E13 .* E24));

varE13 = E13sq_m - E13m^2;
varE24 = E24sq_m - E24m^2;
covE13E24 = E13E24_m - E13m * E24m;

de13_dT13 = varE13 / (AD.k^2 * t13^2);
de13_dT24 = covE13E24 / (AD.k^2 * t24^2);
de24_dT13 = covE13E24 / (AD.k^2 * t13^2);
de24_dT24 = varE24 / (AD.k^2 * t24^2);

A = zeros(3, 3);
A(1, 1) = 3;
A(1, 2) = de13_dT13 + de24_dT13;
A(1, 3) = de13_dT24 + de24_dT24;
A(2, 2) = de13_dT13;
A(2, 3) = de13_dT24;
A(3, 2) = de24_dT13;
A(3, 3) = de24_dT24;

B = [0; R13; R24];
dy = A \ B;

end

function kb = k_vt_bwd(kf, i, f, T, AD)
kb = kf * (AD.stw(i) / AD.stw(f)) * exp((AD.e1234(f) - AD.e1234(i)) / (AD.k * T));
end

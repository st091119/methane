function dy = rpart_mt_3t_sts(~, y, AD)
% RPART_MT_3T_STS  Трехтемпературная CH4-модель: y = [T_b; T13_b; T24_b].
% Источники считаются по state-to-state каналам VT2, VT4, VV34, VV34d, VV34_4d.

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

if ~isfield(AD, 'fho_steric_vt_v3')
    AD.fho_steric_vt_v3 = sqrt(AD.fho_steric2 * AD.fho_steric4);
end
if ~isfield(AD, 'fho_steric_vv_34s')
    AD.fho_steric_vv_34s = 0.068;
end
if ~isfield(AD, 'fho_steric_vt_v3_4')
    AD.fho_steric_vt_v3_4 = 0.99;
end
if ~isfield(AD, 'fho_steric_vv_34d')
    AD.fho_steric_vv_34d = 1.0;
end
if ~isfield(AD, 'fho_steric_vv_34_4d')
    AD.fho_steric_vv_34_4d = 0.2;
end

n_scale = AD.n0 * AD.tau;
sk0 = [0, 0, 0, 0];

Jmx = AD.lch4(2) - 1;
Lmx = AD.lch4(4) - 1;

k2_tab = zeros(Jmx + 1, 1);
for Jq = 0:Jmx - 1
    k2_tab(Jq + 1) = rates_fho_vt(temp, Jq, Jq + 1, 2, AD.fho_steric2, AD.fho_alpha2, AD.fho_e_m2, 'n_steps', 8000);
end

k4_tab = zeros(Lmx + 1, 1);
for Lq = 0:Lmx - 1
    k4_tab(Lq + 1) = rates_fho_vt(temp, Lq, Lq + 1, 4, AD.fho_steric4, AD.fho_alpha4, AD.fho_e_m4, 'n_steps', 8000);
end

RVT2 = zeros(N, 1);
src = find(~isnan(AD.indj_vt2));
for ii = 1:numel(src)
    r = src(ii);
    dst = AD.indj_vt2(r);
    Jq = AD.inds(r, 2);
    kf = k2_tab(Jq + 1) * n_scale;
    kr = k_vt_bwd(kf, r, dst, temp, AD);
    RVT2(r) = RVT2(r) + nco2i_b(dst) * kr - nco2i_b(r) * kf;
    RVT2(dst) = RVT2(dst) + nco2i_b(r) * kf - nco2i_b(dst) * kr;
end

RVT4 = zeros(N, 1);
src = find(~isnan(AD.indl_vt4));
for ii = 1:numel(src)
    r = src(ii);
    dst = AD.indl_vt4(r);
    Lq = AD.inds(r, 4);
    kf = k4_tab(Lq + 1) * n_scale;
    kr = k_vt_bwd(kf, r, dst, temp, AD);
    RVT4(r) = RVT4(r) + nco2i_b(dst) * kr - nco2i_b(r) * kf;
    RVT4(dst) = RVT4(dst) + nco2i_b(r) * kf - nco2i_b(dst) * kr;
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
        kv = rates_fho_vv(temp, si, sf, sk0, sk0, AD.fho_steric_vt_v3, AD.fho_steric_vv_34s, AD.fho_alpha4, AD.fho_e_m4, 'n_steps', 6000);
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
        kv = rates_fho_vv(temp, si, sf, sk0, sk0, AD.fho_steric_vt_v3_4, AD.fho_steric_vv_34d, AD.fho_alpha4, AD.fho_e_m4, 'n_steps', 6000);
        kvv34d_cache(key) = kv;
    end
    kf = kv * n_scale;
    kr = k_vt_bwd(kf, r, dst, temp, AD);
    RVV34d(r) = RVV34d(r) + nco2i_b(dst) * kr - nco2i_b(r) * kf;
    RVV34d(dst) = RVV34d(dst) + nco2i_b(r) * kf - nco2i_b(dst) * kr;
end

RVV34_4d = zeros(N, 1);
src = find(~isnan(AD.indvv34));
kvv34_4d_cache = containers.Map('KeyType', 'char', 'ValueType', 'double');
for ii = 1:numel(src)
    r = src(ii);
    dst = AD.indvv34(r);
    si = AD.inds(r, :);
    sf = AD.inds(dst, :);
    key = sprintf('%d_%d', si(3), si(4));
    if isKey(kvv34_4d_cache, key)
        kv = kvv34_4d_cache(key);
    else
        kv = rates_fho_vv(temp, si, sf, sk0, [0, 0, 0, 1], AD.fho_steric_vt_v3, AD.fho_steric_vv_34_4d, AD.fho_alpha4, AD.fho_e_m4, 'n_steps', 6000);
        kvv34_4d_cache(key) = kv;
    end
    kf = kv * n_scale;
    kr = k_vt_bwd(kf, r, dst, temp, AD);
    RVV34_4d(r) = RVV34_4d(r) + nco2i_b(dst) * kr - nco2i_b(r) * kf;
    RVV34_4d(dst) = RVV34_4d(dst) + nco2i_b(r) * kf - nco2i_b(dst) * kr;
end

Rtot = RVT2 + RVT4 + RVV + RVV34d + RVV34_4d;
R13 = sum((E13 / kT0) .* Rtot);
R24 = sum((E24 / kT0) .* Rtot);

% Численные производные средних энергий по T13 и T24 для матрицы A.
h = 1e-4;
[e13_p, e24_p] = grouped_energies(T13_b + h, T24_b, AD, E13, E24);
[e13_m, e24_m] = grouped_energies(max(T13_b - h, 1e-9), T24_b, AD, E13, E24);
de13_dT13 = (e13_p - e13_m) / ((T13_b + h) - max(T13_b - h, 1e-9));
de24_dT13 = (e24_p - e24_m) / ((T13_b + h) - max(T13_b - h, 1e-9));

[e13_p, e24_p] = grouped_energies(T13_b, T24_b + h, AD, E13, E24);
[e13_m, e24_m] = grouped_energies(T13_b, max(T24_b - h, 1e-9), AD, E13, E24);
de13_dT24 = (e13_p - e13_m) / ((T24_b + h) - max(T24_b - h, 1e-9));
de24_dT24 = (e24_p - e24_m) / ((T24_b + h) - max(T24_b - h, 1e-9));

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

function [e13_nd, e24_nd] = grouped_energies(T13_b, T24_b, AD, E13, E24)
t13 = max(T13_b, 1e-9) * AD.T0;
t24 = max(T24_b, 1e-9) * AD.T0;
xi = -(E13 ./ (AD.k * t13) + E24 ./ (AD.k * t24));
f = AD.stw(:) .* exp(xi);
Zv = sum(f);
nbi = f ./ Zv;
kT0 = AD.k * AD.T0;
e13_nd = sum(nbi .* E13) / kT0;
e24_nd = sum(nbi .* E24) / kT0;
end

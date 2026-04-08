function dy = rpart_mt_sts(~, y, AD)

%   y = [T_b; Tv_b], где T_b = T/T0, Tv_b = Tv/T0.
%   Процессы: VT2, VT4, VV34s, VV34d, VV34_4d.

if ~exist('rates_fho_vt', 'file')
    addpath(fullfile(fileparts(mfilename('fullpath')), 'fho_model'));
end

AD = ch4_relax_topology(AD);

%% безразмерные температуры (нижняя граница для устойчивости)
T_b  = max(y(1), 1e-9);
Tv_b = max(y(2), 1e-9);

num_var = 2;
kT0 = AD.k * AD.T0;

temp = T_b * AD.T0;
tv = Tv_b * AD.T0;

%% колебательное распределение (как в vdf.m / rpart_mt_lt.m)
w = AD.stw(:);
E = AD.e1234(:);
fac = -E ./ (AD.k * tv);
f = w .* exp(fac);
Zv = sum(f);
nco2i_b = f ./ Zv;

%% стерические факторы VV (по умолчанию как в fho_model/test_fho.m)
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

sk0 = [0, 0, 0, 0];
n_scale = AD.n0 * AD.tau;

N = numel(AD.e1234);
qnco2 = N;

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

%% --- VT2 ---
RVT2 = zeros(qnco2, 1);
indj = AD.indj_vt2;
vt2_src = find(~isnan(indj));

for ii = 1:numel(vt2_src)
    r = vt2_src(ii);
    dst = indj(r); % (I,J,K,L) -> (I,J+1,K,L)
    Jq = AD.inds(r, 2);
    kf = k2_tab(Jq + 1) * n_scale;
    kr = k_vt_bwd(kf, r, dst, temp, AD);
    RVT2(r) = RVT2(r) + nco2i_b(dst) * kr - nco2i_b(r) * kf;
    RVT2(dst) = RVT2(dst) + nco2i_b(r) * kf - nco2i_b(dst) * kr;
end

RVIBR = sum((E / kT0) .* RVT2);

%% --- VT4 ---
RVT4 = zeros(qnco2, 1);
indl = AD.indl_vt4;
vt4_src = find(~isnan(indl));

for ii = 1:numel(vt4_src)
    r = vt4_src(ii);
    dst = indl(r); % (I,J,K,L) -> (I,J,K,L+1)
    Lq = AD.inds(r, 4);
    kf = k4_tab(Lq + 1) * n_scale;
    kr = k_vt_bwd(kf, r, dst, temp, AD);
    RVT4(r) = RVT4(r) + nco2i_b(dst) * kr - nco2i_b(r) * kf;
    RVT4(dst) = RVT4(dst) + nco2i_b(r) * kf - nco2i_b(dst) * kr;
end

RVIBR = RVIBR + sum((E / kT0) .* RVT4);

%% --- VV34: (K,L) <-> (K-1, L+1), внутримолекулярный FHO ---
% Обратный переход (K-1,L+1)->(K,L), т.е. "K+1,L-1" относительно dst,
% учитывается через k_r = detailed balance(k_f).
RVV = zeros(qnco2, 1);
indvv34 = AD.indvv34;
vv_src = find(~isnan(indvv34));
kvv_up = zeros(N, 1);
% Кэш FHO по (K,L) донора: k(I,J,K,L) слабо зависит от I,J относительно полного перебора уровней.
kvv_cache = containers.Map('KeyType', 'char', 'ValueType', 'double');
for ii = 1:numel(vv_src)
    r = vv_src(ii);
    dst = indvv34(r);
    si = AD.inds(r, :);
    sf = AD.inds(dst, :);
    key = sprintf('%d_%d', si(3), si(4));
    if isKey(kvv_cache, key)
        kvv_up(r) = kvv_cache(key);
    else
        kv = rates_fho_vv(temp, si, sf, sk0, sk0, AD.fho_steric_vt_v3, AD.fho_steric_vv_34s, ...
            AD.fho_alpha4, AD.fho_e_m4, 'n_steps', 6000);
        kvv_cache(key) = kv;
        kvv_up(r) = kv;
    end
end
kvv_up = kvv_up * n_scale;

for ii = 1:numel(vv_src)
    r = vv_src(ii);
    dst = indvv34(r);
    kf = kvv_up(r);
    kr = k_vt_bwd(kf, r, dst, temp, AD);
    RVV(r) = RVV(r) + nco2i_b(dst) * kr - nco2i_b(r) * kf;
    RVV(dst) = RVV(dst) + nco2i_b(r) * kf - nco2i_b(dst) * kr;
end

RVIBR = RVIBR + sum((E / kT0) .* RVV);

%% --- VV34d: CH4(i) + CH4 <-> CH4(i3-1,i4+2) + CH4 ---
% Обратный канал i3+1,i4-2 также учитывается через k_r (detailed balance).
RVV34d = zeros(qnco2, 1);
indvv34d = AD.indvv34d;
vv34d_src = find(~isnan(indvv34d));
kvv34d_up = zeros(N, 1);
kvv34d_cache = containers.Map('KeyType', 'char', 'ValueType', 'double');
for ii = 1:numel(vv34d_src)
    r = vv34d_src(ii);
    dst = indvv34d(r);
    si = AD.inds(r, :);
    sf = AD.inds(dst, :);
    key = sprintf('%d_%d', si(3), si(4));
    if isKey(kvv34d_cache, key)
        kvv34d_up(r) = kvv34d_cache(key);
    else
        kv = rates_fho_vv(temp, si, sf, sk0, sk0, AD.fho_steric_vt_v3_4, AD.fho_steric_vv_34d, ...
            AD.fho_alpha4, AD.fho_e_m4, 'n_steps', 6000);
        kvv34d_cache(key) = kv;
        kvv34d_up(r) = kv;
    end
end
kvv34d_up = kvv34d_up * n_scale;
for ii = 1:numel(vv34d_src)
    r = vv34d_src(ii);
    dst = indvv34d(r);
    kf = kvv34d_up(r);
    kr = k_vt_bwd(kf, r, dst, temp, AD);
    RVV34d(r) = RVV34d(r) + nco2i_b(dst) * kr - nco2i_b(r) * kf;
    RVV34d(dst) = RVV34d(dst) + nco2i_b(r) * kf - nco2i_b(dst) * kr;
end
RVIBR = RVIBR + sum((E / kT0) .* RVV34d);

%% --- VV34_4d: CH4(i)+CH4(k) <-> CH4(i3-1,i4+1)+CH4(k4+1) ---
% Обратный канал CH4(i3+1,i4-1)+CH4(k4-1) идет через k_r.
RVV34_4d = zeros(qnco2, 1);
vv34_4d_src = find(~isnan(indvv34));
kvv34_4d_up = zeros(N, 1);
kvv34_4d_cache = containers.Map('KeyType', 'char', 'ValueType', 'double');
for ii = 1:numel(vv34_4d_src)
    r = vv34_4d_src(ii);
    dst = indvv34(r);
    si = AD.inds(r, :);
    sf = AD.inds(dst, :);
    key = sprintf('%d_%d', si(3), si(4));
    if isKey(kvv34_4d_cache, key)
        kvv34_4d_up(r) = kvv34_4d_cache(key);
    else
        kv = rates_fho_vv(temp, si, sf, sk0, [0, 0, 0, 1], AD.fho_steric_vt_v3, AD.fho_steric_vv_34_4d, ...
            AD.fho_alpha4, AD.fho_e_m4, 'n_steps', 6000);
        kvv34_4d_cache(key) = kv;
        kvv34_4d_up(r) = kv;
    end
end
kvv34_4d_up = kvv34_4d_up * n_scale;
for ii = 1:numel(vv34_4d_src)
    r = vv34_4d_src(ii);
    dst = indvv34(r);
    kf = kvv34_4d_up(r);
    kr = k_vt_bwd(kf, r, dst, temp, AD);
    RVV34_4d(r) = RVV34_4d(r) + nco2i_b(dst) * kr - nco2i_b(r) * kf;
    RVV34_4d(dst) = RVV34_4d(dst) + nco2i_b(r) * kf - nco2i_b(dst) * kr;
end
RVIBR = RVIBR + sum((E / kT0) .* RVV34_4d);

%% матрица A для двухтемпературной постановки
e_sum_e_sum = sum(w .* (E ./ (AD.k * tv)) .* exp(fac))^2;
e_sum_square = sum(w .* (E ./ (AD.k * tv)).^2 .* exp(fac));

A = eye(num_var);
A(1, 1) = 3;
A(1, 2) = e_sum_square / Zv - e_sum_e_sum / (Zv^2);
A(2, 2) = A(1, 2);

B = zeros(num_var, 1);
B(2) = RVIBR;

dy = A \ B;

end

% --- detailed balance: обратный к i->f при заданном k_{i->f} ---
function kb = k_vt_bwd(kf, i, f, T, AD)
kb = kf * (AD.stw(i) / AD.stw(f)) * exp((AD.e1234(f) - AD.e1234(i)) / (AD.k * T));
end

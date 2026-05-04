function dy = rpart_mt_sts(~, y, AD)

%   y = [T_b; Tv_b], где T_b = T/T0, Tv_b = Tv/T0.
%   Процессы: VT2, VT4, VV34s, VV34d, VV32d (VV^d_{3-2}), VV34_4d.

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
if ~isfield(AD, 'fho_steric_vv_32d')
    AD.fho_steric_vv_32d = 0.0025;
end
if ~isfield(AD, 'fho_gamma_vv32d')
    AD.fho_gamma_vv32d = 0.5;
end
if ~isfield(AD, 'sw_vv34d_model')
    AD.sw_vv34d_model = 'hard'; % 'hard' | 'easy' (legacy aliases: 'state' | 'macro')
end
if ~isfield(AD, 'sw_vv34_4d_model')
    AD.sw_vv34_4d_model = 'hard'; % 'hard' | 'easy' (legacy aliases: 'state' | 'macro')
end
if ~isfield(AD, 'n_steps_vv_easy')
    AD.n_steps_vv_easy = 300;
end
if ~isfield(AD, 'use_vv34')
    AD.use_vv34 = true;
end
if ~isfield(AD, 'use_vv34d')
    AD.use_vv34d = true;
end
if ~isfield(AD, 'use_vv32d')
    AD.use_vv32d = true;
end
if ~isfield(AD, 'use_vv34_4d')
    AD.use_vv34_4d = true;
end

sk0 = [0, 0, 0, 0];
n_scale = AD.n0 * AD.tau;

N = numel(AD.e1234);
qnco2 = N;

Jmx = AD.lch4(2) - 1;
Lmx = AD.lch4(4) - 1;

% VT2: deexcitation rates k_{J -> J-1} stored at index J (for J = 1..Jmx)
k2_tab = zeros(Jmx + 1, 1);
for Jq = 1:Jmx
    k2_tab(Jq) = rates_fho_vt(temp, Jq, Jq - 1, 2, AD.fho_steric2, AD.fho_alpha, AD.fho_e_m, 'n_steps', 8000);
end

% VT4: deexcitation rates k_{L -> L-1} stored at index L (for L = 1..Lmx)
k4_tab = zeros(Lmx + 1, 1);
for Lq = 1:Lmx
    k4_tab(Lq) = rates_fho_vt(temp, Lq, Lq - 1, 4, AD.fho_steric4, AD.fho_alpha, AD.fho_e_m, 'n_steps', 8000);
end

%% --- VT2 ---
RVT2 = zeros(qnco2, 1);
indj = AD.indj_vt2;
vt2_src = find(~isnan(indj));

for ii = 1:numel(vt2_src)
    r = vt2_src(ii);
    dst = indj(r); % low J -> high J neighbor
    Jq = AD.inds(r, 2);
    kf = k2_tab(Jq + 1) * n_scale;      % forward: dst(J+1) -> r(J)
    kr = k_vt_bwd(kf, dst, r, temp, AD); % backward: r(J) -> dst(J+1)
    RVT2(r) = RVT2(r) + nco2i_b(dst) * kf - nco2i_b(r) * kr;
    RVT2(dst) = RVT2(dst) + nco2i_b(r) * kr - nco2i_b(dst) * kf;
end

RVIBR = sum((E / kT0) .* RVT2);

%% --- VT4 ---
RVT4 = zeros(qnco2, 1);
indl = AD.indl_vt4;
vt4_src = find(~isnan(indl));

for ii = 1:numel(vt4_src)
    r = vt4_src(ii);
    dst = indl(r); % low L -> high L neighbor
    Lq = AD.inds(r, 4);
    kf = k4_tab(Lq + 1) * n_scale;      % forward: dst(L+1) -> r(L)
    kr = k_vt_bwd(kf, dst, r, temp, AD); % backward: r(L) -> dst(L+1)
    RVT4(r) = RVT4(r) + nco2i_b(dst) * kf - nco2i_b(r) * kr;
    RVT4(dst) = RVT4(dst) + nco2i_b(r) * kr - nco2i_b(dst) * kf;
end

RVIBR = RVIBR + sum((E / kT0) .* RVT4);

%% --- VV34: (K,L) <-> (K-1, L+1), внутримолекулярный FHO ---
% Обратный переход (K-1,L+1)->(K,L), т.е. "K+1,L-1" относительно dst,
% учитывается через k_r = detailed balance(k_f).
RVV = zeros(qnco2, 1);
if AD.use_vv34
indvv34 = AD.indvv34;
vv_src = find(~isnan(indvv34));
kvv_up = zeros(N, 1);
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
        kv = rates_fho_vv(temp, si, sf, sk0, sk0, AD.fho_steric_vv_34s, ...
            AD.fho_alpha, AD.fho_e_m, 'n_steps', 6000);
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
end

RVIBR = RVIBR + sum((E / kT0) .* RVV);

%% --- VV34d: CH4(i) + CH4 <-> CH4(i3-1,i4+2) + CH4 ---
if AD.use_vv34d
vv34d_is_easy = strcmpi(AD.sw_vv34d_model, 'easy') || strcmpi(AD.sw_vv34d_model, 'macro');
if vv34d_is_easy
    n_steps_vv34d = AD.n_steps_vv_easy;
    A34d = ch4_A_vv34d(temp, tv, tv, AD, n_steps_vv34d);
    eps3_vv = AD.e0010 - AD.e0000;
    eps4_vv = AD.e0001 - AD.e0000;
    RVIBR = RVIBR + (AD.n0 * A34d * n_scale / kT0) * (2 * eps4_vv - eps3_vv);
else
    % State-resolved legacy branch.
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
            kv = rates_fho_vv(temp, si, sf, sk0, sk0, AD.fho_steric_vv_34d, ...
                AD.fho_alpha, AD.fho_e_m, 'n_steps', 6000, ...
                'gamma', AD.fho_gamma_vv34d);
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
end
end

%% --- VV32d: CH4(i) + CH4 <-> CH4(i2+2,i3-1) + CH4, VV^d_{3-2} ---
RVV32d = zeros(qnco2, 1);
if AD.use_vv32d
indvv32d = AD.indvv32d;
vv32d_src = find(~isnan(indvv32d));
kvv32d_up = zeros(N, 1);
kvv32d_cache = containers.Map('KeyType', 'char', 'ValueType', 'double');
for ii = 1:numel(vv32d_src)
    r = vv32d_src(ii);
    dst = indvv32d(r);
    si = AD.inds(r, :);
    sf = AD.inds(dst, :);
    key = sprintf('%d_%d', si(2), si(3));
    if isKey(kvv32d_cache, key)
        kvv32d_up(r) = kvv32d_cache(key);
    else
        kv = rates_fho_vv(temp, si, sf, sk0, sk0, AD.fho_steric_vv_32d, ...
            AD.fho_alpha, AD.fho_e_m, 'n_steps', 6000, ...
            'gamma', AD.fho_gamma_vv32d);
        kvv32d_cache(key) = kv;
        kvv32d_up(r) = kv;
    end
end
kvv32d_up = kvv32d_up * n_scale;
for ii = 1:numel(vv32d_src)
    r = vv32d_src(ii);
    dst = indvv32d(r);
    kf = kvv32d_up(r);
    kr = k_vt_bwd(kf, r, dst, temp, AD);
    RVV32d(r) = RVV32d(r) + nco2i_b(dst) * kr - nco2i_b(r) * kf;
    RVV32d(dst) = RVV32d(dst) + nco2i_b(r) * kf - nco2i_b(dst) * kr;
end
end
RVIBR = RVIBR + sum((E / kT0) .* RVV32d);

%% --- VV34_4d: CH4(i)+CH4(k4) <-> CH4(i3-1,i4+1)+CH4(k4+1) ---
if AD.use_vv34_4d
vv34_4d_is_easy = strcmpi(AD.sw_vv34_4d_model, 'easy') || strcmpi(AD.sw_vv34_4d_model, 'macro');
if vv34_4d_is_easy
    n_steps_vv34_4d = AD.n_steps_vv_easy;
    B34 = ch4_B_vv34_4d(temp, tv, tv, AD, n_steps_vv34_4d);
    eps3_vv = AD.e0010 - AD.e0000;
    eps4_vv = AD.e0001 - AD.e0000;
    RVIBR = RVIBR + (AD.n0 * B34 * n_scale / kT0) * (2 * eps4_vv - eps3_vv);
else
    [RVV34_4d, RVV34_4d_partner_E] = ch4_vv34_4d_source(temp, nco2i_b, AD, n_scale, 6000);
    RVIBR = RVIBR + sum((E / kT0) .* RVV34_4d) + RVV34_4d_partner_E / kT0;
end
end

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

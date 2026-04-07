function dy = rpart_mt_sts(~, y, AD)
%RPART_MT_STS  Правые части СОДУ: поуровневые VT2, VT4 (FHO) и VV34 (ν₃↔ν₄, 1 квант).
%   y = [T_b; Tv2_b; Tv3_b; Tv4_b] — безразмерные температуры (масштаб AD.T0).
%   Мода ν₁ описывается поступательной T; моды 2–4 — своими Tv2,Tv3,Tv4.
%
%   Требуются в AD: поля из input_data + n0,T0,tau + параметры FHO (как в main.m):
%   fho_alpha2, fho_e_m2, fho_alpha4, fho_e_m4, fho_steric2, fho_steric4.
%   Опционально: fho_steric_vt_v3 (по умолч. sqrt(steric2*steric4)),
%              fho_steric_vv_34s (по умолч. 0.068, как VV34s в test_fho.m).

if ~exist('rates_fho_vt', 'file')
    addpath(fullfile(fileparts(mfilename('fullpath')), 'fho_model'));
end

AD = ch4_relax_topology(AD);

%% безразмерные температуры (нижняя граница для устойчивости)
T_b   = max(y(1), 1e-9);
Tv2_b = max(y(2), 1e-9);
Tv3_b = max(y(3), 1e-9);
Tv4_b = max(y(4), 1e-9);

num_var = 4;
kT0 = AD.k * AD.T0;

temp = T_b * AD.T0;
t2 = Tv2_b * AD.T0;
t3 = Tv3_b * AD.T0;
t4 = Tv4_b * AD.T0;

%% частичные энергии [Дж] и безразмерные e / (k T0)
hc = AD.h * AD.c;
d = AD.d;
I0 = AD.inds(:, 1);
J0 = AD.inds(:, 2);
K0 = AD.inds(:, 3);
L0 = AD.inds(:, 4);

E1p = hc .* AD.omega(1) .* (I0 + d(1)/2);
E2p = hc .* AD.omega(2) .* (J0 + d(2)/2);
E3p = hc .* AD.omega(3) .* (K0 + d(3)/2);
E4p = hc .* AD.omega(4) .* (L0 + d(4)/2);

e1_b = E1p / kT0;
e2_b = E2p / kT0;
e3_b = E3p / kT0;
e4_b = E4p / kT0;
Ev_b = AD.e1234 / kT0;

xi = -e1_b ./ T_b - e2_b ./ Tv2_b - e3_b ./ Tv3_b - e4_b ./ Tv4_b;
f = AD.stw .* exp(xi);
Zv = sum(f);
nco2i_b = f ./ Zv;

%% стерические факторы VV (по умолчанию как в fho_model/test_fho.m)
if ~isfield(AD, 'fho_steric_vt_v3')
    AD.fho_steric_vt_v3 = sqrt(AD.fho_steric2 * AD.fho_steric4);
end
if ~isfield(AD, 'fho_steric_vv_34s')
    AD.fho_steric_vv_34s = 0.068;
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
indjp = AD.indjp_vt2;
indjm = AD.indjm_vt2;
jp = find(~isnan(indjp));
jm = find(~isnan(indjm));

kvt2_up = zeros(N, 1);
for ii = 1:numel(jp)
    r = jp(ii);
    Jq = AD.inds(r, 2);
    kvt2_up(r) = k2_tab(Jq + 1);
end
kvt2_up = kvt2_up * n_scale;

for ii = 1:numel(jp)
    r = jp(ii);
    r1 = indjp(r);
    kvt2_down_r = k_vt_bwd(kvt2_up(r), r, r1, temp, AD);
    RVT2(r) = RVT2(r) + nco2i_b(r1) * kvt2_down_r - nco2i_b(r) * kvt2_up(r);
end
for ii = 1:numel(jm)
    r = jm(ii);
    r1 = indjm(r);
    k_up_r1 = kvt2_up(r1);
    if k_up_r1 == 0
        continue
    end
    k_dn_r = k_vt_bwd(k_up_r1, r1, r, temp, AD);
    RVT2(r) = RVT2(r) + nco2i_b(r1) * k_up_r1 - nco2i_b(r) * k_dn_r;
end

RVIB_2 = sum(E2p / kT0 .* RVT2);

%% --- VT4 ---
RVT4 = zeros(qnco2, 1);
jp4 = find(~isnan(AD.indjp_vt4));
jm4 = find(~isnan(AD.indjm_vt4));

kvt4_up = zeros(N, 1);
for ii = 1:numel(jp4)
    r = jp4(ii);
    Lq = AD.inds(r, 4);
    kvt4_up(r) = k4_tab(Lq + 1);
end
kvt4_up = kvt4_up * n_scale;

for ii = 1:numel(jp4)
    r = jp4(ii);
    r1 = AD.indjp_vt4(r);
    kdn = k_vt_bwd(kvt4_up(r), r, r1, temp, AD);
    RVT4(r) = RVT4(r) + nco2i_b(r1) * kdn - nco2i_b(r) * kvt4_up(r);
end
for ii = 1:numel(jm4)
    r = jm4(ii);
    r1 = AD.indjm_vt4(r);
    k_up_r1 = kvt4_up(r1);
    if k_up_r1 == 0
        continue
    end
    k_dn_r = k_vt_bwd(k_up_r1, r1, r, temp, AD);
    RVT4(r) = RVT4(r) + nco2i_b(r1) * k_up_r1 - nco2i_b(r) * k_dn_r;
end

RVIB_4 = sum(E4p / kT0 .* RVT4);

%% --- VV34: (K,L) -> (K-1, L+1), внутримолекулярный FHO (partner sk0 неизменен) ---
RVV = zeros(qnco2, 1);
vv_src = find(~isnan(AD.indvv34_dst));
kvv_up = zeros(N, 1);
% Кэш FHO по (K,L) донора: k(I,J,K,L) слабо зависит от I,J относительно полного перебора уровней.
kvv_cache = containers.Map('KeyType', 'char', 'ValueType', 'double');
for ii = 1:numel(vv_src)
    r = vv_src(ii);
    dst = AD.indvv34_dst(r);
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
    dst = AD.indvv34_dst(r);
    kf = kvv_up(r);
    kr = k_vt_bwd(kf, r, dst, temp, AD);
    RVV(r) = RVV(r) + nco2i_b(dst) * kr - nco2i_b(r) * kf;
    RVV(dst) = RVV(dst) + nco2i_b(r) * kf - nco2i_b(dst) * kr;
end

RVIB_3 = sum(E3p / kT0 .* RVV);
RVIB_4 = RVIB_4 + sum(E4p / kT0 .* RVV);

%% матрица A: строка 1 — баланс энергии (3 — поступательная + вращательная на моль); 2–4 — средние энергии мод
A = zeros(num_var, num_var);
row_ev = dmu_row(Ev_b, f, Zv, e1_b, e2_b, e3_b, e4_b, T_b, Tv2_b, Tv3_b, Tv4_b);
row_e2 = dmu_row(E2p / kT0, f, Zv, e1_b, e2_b, e3_b, e4_b, T_b, Tv2_b, Tv3_b, Tv4_b);
row_e3 = dmu_row(E3p / kT0, f, Zv, e1_b, e2_b, e3_b, e4_b, T_b, Tv2_b, Tv3_b, Tv4_b);
row_e4 = dmu_row(E4p / kT0, f, Zv, e1_b, e2_b, e3_b, e4_b, T_b, Tv2_b, Tv3_b, Tv4_b);

A(1, :) = [3, 0, 0, 0] + row_ev;
A(2, :) = row_e2;
A(3, :) = row_e3;
A(4, :) = row_e4;

B = zeros(num_var, 1);
B(2) = RVIB_2;
B(3) = RVIB_3;
B(4) = RVIB_4;

dy = A \ B;

end

% --- detailed balance: обратный к i->f при заданном k_{i->f} ---
function kb = k_vt_bwd(kf, i, f, T, AD)
kb = kf * (AD.stw(i) / AD.stw(f)) * exp((AD.e1234(f) - AD.e1234(i)) / (AD.k * T));
end

function row = dmu_row(c, f, Zv, e1_b, e2_b, e3_b, e4_b, Tb, T2, T3, T4)
dfdTb = f .* (e1_b ./ (Tb .^ 2));
dfdT2 = f .* (e2_b ./ (T2 .^ 2));
dfdT3 = f .* (e3_b ./ (T3 .^ 2));
dfdT4 = f .* (e4_b ./ (T4 .^ 2));
dZdTb = sum(dfdTb);
dZdT2 = sum(dfdT2);
dZdT3 = sum(dfdT3);
dZdT4 = sum(dfdT4);
num = sum(c .* f);
row = zeros(1, 4);
row(1) = (sum(c .* dfdTb) * Zv - num * dZdTb) / (Zv ^ 2);
row(2) = (sum(c .* dfdT2) * Zv - num * dZdT2) / (Zv ^ 2);
row(3) = (sum(c .* dfdT3) * Zv - num * dZdT3) / (Zv ^ 2);
row(4) = (sum(c .* dfdT4) * Zv - num * dZdT4) / (Zv ^ 2);
end

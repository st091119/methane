clear

% Analysis of CH4 process contributions Q^gamma [J/(m^3 sec)].
% Each process gamma is evaluated along the 3T STS relaxation trajectory.

addpath(fullfile(fileparts(mfilename('fullpath')), 'fho_model'));

%% Input parameters
p0 = 101325;        % pressure [Pa]
T0 = 1000;           % translational temperature [K]
T13_0 = 700;        % temperature for modes 1 and 3 [K]
T24_0 = 300;        % temperature for modes 2 and 4 [K]
t_fin = 1e-2;       % final time for contribution plot [sec]

options = odeset('RelTol', 1e-8, 'AbsTol', 1e-8);

%% Model data
AD = input_data;

n0 = p0 / (AD.k * T0);
sigma0 = pi * AD.r0^2;
tau = (4 * n0 * sigma0 * sqrt(AD.k * T0 / (pi * AD.m)))^(-1);

AD.n0 = n0;
AD.T0 = T0;
AD.p0 = p0;
AD.tau = tau;

AD.fho_alpha   = 5.174e10;
AD.fho_e_m     = 776.4;
AD.fho_steric2 = 0.002612;
AD.fho_steric4 = 0.028546;

AD.fho_steric_vv_34s = 1.0;
AD.fho_steric_vv_34d = 0.24;
AD.fho_steric_vv_32d = 0.0025;
AD.fho_gamma_vv32d = 0.5;
AD.fho_gamma_vv34d = 1.0;
AD.fho_steric_vv_34_4d = 1.0;
AD.fho_steric_vv_32_2d = 1.0;
AD.fho_steric_vv_34_2d = 1.0;
AD.fho_steric_vv_32_4d = 1.0;

AD.sw_vv34d_model = 'hard';
AD.sw_vv34_4d_model = 'hard';
AD.n_steps_vv_easy = 1000;

AD.use_vv34 = true;
AD.use_vv34d = true;
AD.use_vv32d = true;
AD.use_vv34_4d = true;
AD.use_vv32_2d = true;
AD.use_vv34_2d = true;
AD.use_vv32_4d = true;

AD = ch4_relax_topology(AD);

%% Solve 3T relaxation
time = logspace(-8, log10(t_fin), 90).';
tspan = [0; time] ./ tau;
Y0 = [1; T13_0 / T0; T24_0 / T0];

[X, Y] = ode15s(@(t, y) rpart_mt_3t_sts(t, y, AD), tspan, Y0, options);

time = X(2:end) * tau;
Y = Y(2:end, :);

%% Compute Q^gamma for each process
process_names = {'VT2', 'VT4', 'VV3-4', 'VV3-4d', 'VV3-2d', ...
    'VV3-4-4d', 'VV3-2-2d', 'VV3-4-2d', 'VV3-2-4d'};
Q = zeros(numel(time), numel(process_names));

for it = 1:numel(time)
    T = Y(it, 1) * T0;
    T13 = Y(it, 2) * T0;
    T24 = Y(it, 3) * T0;
    contrib = ch4_process_Q_3t(T, T13, T24, AD);
    Q(it, :) = [contrib.VT2, contrib.VT4, contrib.VV34, contrib.VV34d, ...
        contrib.VV32d, contrib.VV34_4d, contrib.VV32_2d, ...
        contrib.VV34_2d, contrib.VV32_4d];
end

%% Plot
colors = [
    1.00, 0.00, 0.00
    0.00, 0.25, 0.70
    1.00, 0.65, 0.00
    0.00, 0.65, 0.25
    0.35, 0.30, 0.60
    0.50, 0.50, 0.50
    0.85, 0.20, 0.75
    0.10, 0.70, 0.80
    0.45, 0.25, 0.10
];

figure;
hold on;
for ip = 1:numel(process_names)
    semilogy(time, max(abs(Q(:, ip)), realmin), 'LineWidth', 2, ...
        'Color', colors(ip, :));
end
set(gca, 'XScale', 'log', 'YScale', 'log');
xlabel('time [sec]');
ylabel('abs(Q^\gamma) [J/(m^3 sec)]');
legend(process_names, 'Location', 'best');
title('CH4 process contributions');
grid on;
box on;

%% Local functions
function Q = ch4_process_Q_3t(T, T13, T24, AD)
if ~exist('rates_fho_vt', 'file')
    addpath(fullfile(fileparts(mfilename('fullpath')), 'fho_model'));
end

T = max(T, 1e-9);
T13 = max(T13, 1e-9);
T24 = max(T24, 1e-9);

kT0 = AD.k * AD.T0;
n_scale = AD.n0 * AD.tau;
dim_scale = AD.n0 * AD.k * AD.T0 / AD.tau;

N = numel(AD.e1234);
sk0 = [0, 0, 0, 0];

eps1 = AD.e1000 - AD.e0000;
eps2 = AD.e0100 - AD.e0000;
eps3 = AD.e0010 - AD.e0000;
eps4 = AD.e0001 - AD.e0000;

I0 = AD.inds(:, 1);
J0 = AD.inds(:, 2);
K0 = AD.inds(:, 3);
L0 = AD.inds(:, 4);
E13 = I0 * eps1 + K0 * eps3;
E24 = J0 * eps2 + L0 * eps4;
E = E13 + E24;

xi = -(E13 ./ (AD.k * T13) + E24 ./ (AD.k * T24));
f = AD.stw(:) .* exp(xi);
x = f ./ sum(f);

Jmx = AD.lch4(2) - 1;
Lmx = AD.lch4(4) - 1;

k2_tab = zeros(Jmx + 1, 1);
for Jq = 1:Jmx
    k2_tab(Jq) = rates_fho_vt(T, Jq, Jq - 1, 2, ...
        AD.fho_steric2, AD.fho_alpha, AD.fho_e_m, 'n_steps', 8000);
end

k4_tab = zeros(Lmx + 1, 1);
for Lq = 1:Lmx
    k4_tab(Lq) = rates_fho_vt(T, Lq, Lq - 1, 4, ...
        AD.fho_steric4, AD.fho_alpha, AD.fho_e_m, 'n_steps', 8000);
end

RVT2 = zeros(N, 1);
src = find(~isnan(AD.indj_vt2));
for ii = 1:numel(src)
    r = src(ii);
    dst = AD.indj_vt2(r);
    Jq = AD.inds(r, 2);
    kf = k2_tab(Jq + 1) * n_scale;
    kr = k_vt_bwd(kf, dst, r, T, AD);
    RVT2(r) = RVT2(r) + x(dst) * kf - x(r) * kr;
    RVT2(dst) = RVT2(dst) + x(r) * kr - x(dst) * kf;
end

RVT4 = zeros(N, 1);
src = find(~isnan(AD.indl_vt4));
for ii = 1:numel(src)
    r = src(ii);
    dst = AD.indl_vt4(r);
    Lq = AD.inds(r, 4);
    kf = k4_tab(Lq + 1) * n_scale;
    kr = k_vt_bwd(kf, dst, r, T, AD);
    RVT4(r) = RVT4(r) + x(dst) * kf - x(r) * kr;
    RVT4(dst) = RVT4(dst) + x(r) * kr - x(dst) * kf;
end

RVV34 = zeros(N, 1);
if AD.use_vv34
    src = find(~isnan(AD.indvv34));
    rate_cache = containers.Map('KeyType', 'char', 'ValueType', 'double');
    for ii = 1:numel(src)
        r = src(ii);
        dst = AD.indvv34(r);
        si = AD.inds(r, :);
        sf = AD.inds(dst, :);
        key = sprintf('%d_%d', si(3), si(4));
        if isKey(rate_cache, key)
            kv = rate_cache(key);
        else
            kv = rates_fho_vv(T, si, sf, sk0, sk0, ...
                AD.fho_steric_vv_34s, AD.fho_alpha, AD.fho_e_m, 'n_steps', 6000);
            rate_cache(key) = kv;
        end
        kf = kv * n_scale;
        kr = k_vt_bwd(kf, r, dst, T, AD);
        RVV34(r) = RVV34(r) + x(dst) * kr - x(r) * kf;
        RVV34(dst) = RVV34(dst) + x(r) * kf - x(dst) * kr;
    end
end

RVV34d = zeros(N, 1);
if AD.use_vv34d
    src = find(~isnan(AD.indvv34d));
    rate_cache = containers.Map('KeyType', 'char', 'ValueType', 'double');
    for ii = 1:numel(src)
        r = src(ii);
        dst = AD.indvv34d(r);
        si = AD.inds(r, :);
        sf = AD.inds(dst, :);
        key = sprintf('%d_%d', si(3), si(4));
        if isKey(rate_cache, key)
            kv = rate_cache(key);
        else
            kv = rates_fho_vv(T, si, sf, sk0, sk0, ...
                AD.fho_steric_vv_34d, AD.fho_alpha, AD.fho_e_m, ...
                'n_steps', 6000, 'gamma', AD.fho_gamma_vv34d);
            rate_cache(key) = kv;
        end
        kf = kv * n_scale;
        kr = k_vt_bwd(kf, r, dst, T, AD);
        RVV34d(r) = RVV34d(r) + x(dst) * kr - x(r) * kf;
        RVV34d(dst) = RVV34d(dst) + x(r) * kf - x(dst) * kr;
    end
end

RVV32d = zeros(N, 1);
if AD.use_vv32d
    src = find(~isnan(AD.indvv32d));
    rate_cache = containers.Map('KeyType', 'char', 'ValueType', 'double');
    for ii = 1:numel(src)
        r = src(ii);
        dst = AD.indvv32d(r);
        si = AD.inds(r, :);
        sf = AD.inds(dst, :);
        key = sprintf('%d_%d', si(2), si(3));
        if isKey(rate_cache, key)
            kv = rate_cache(key);
        else
            kv = rates_fho_vv(T, si, sf, sk0, sk0, ...
                AD.fho_steric_vv_32d, AD.fho_alpha, AD.fho_e_m, ...
                'n_steps', 6000, 'gamma', AD.fho_gamma_vv32d);
            rate_cache(key) = kv;
        end
        kf = kv * n_scale;
        kr = k_vt_bwd(kf, r, dst, T, AD);
        RVV32d(r) = RVV32d(r) + x(dst) * kr - x(r) * kf;
        RVV32d(dst) = RVV32d(dst) + x(r) * kf - x(dst) * kr;
    end
end

if AD.use_vv34_4d
    [RVV34_4d, RVV34_4d_partner_E] = ch4_vv34_4d_source(T, x, AD, n_scale, 6000);
else
    RVV34_4d = zeros(N, 1);
    RVV34_4d_partner_E = 0;
end

if AD.use_vv32_2d
    [RVV32_2d, RVV32_2d_partner_E] = ch4_vv_partner_source(T, x, AD, ...
        n_scale, 6000, [0, 1, -1, 0], 2, AD.fho_steric_vv_32_2d, 0.5);
else
    RVV32_2d = zeros(N, 1);
    RVV32_2d_partner_E = 0;
end

if AD.use_vv34_2d
    [RVV34_2d, RVV34_2d_partner_E] = ch4_vv_partner_source(T, x, AD, ...
        n_scale, 6000, [0, 0, -1, 1], 2, AD.fho_steric_vv_34_2d, 0.5);
else
    RVV34_2d = zeros(N, 1);
    RVV34_2d_partner_E = 0;
end

if AD.use_vv32_4d
    [RVV32_4d, RVV32_4d_partner_E] = ch4_vv_partner_source(T, x, AD, ...
        n_scale, 6000, [0, 1, -1, 0], 4, AD.fho_steric_vv_32_4d, 0.5);
else
    RVV32_4d = zeros(N, 1);
    RVV32_4d_partner_E = 0;
end

Q = struct();
Q.VT2 = dim_scale * sum((E / kT0) .* RVT2);
Q.VT4 = dim_scale * sum((E / kT0) .* RVT4);
Q.VV34 = dim_scale * sum((E / kT0) .* RVV34);
Q.VV34d = dim_scale * sum((E / kT0) .* RVV34d);
Q.VV32d = dim_scale * sum((E / kT0) .* RVV32d);
Q.VV34_4d = dim_scale * (sum((E / kT0) .* RVV34_4d) + RVV34_4d_partner_E / kT0);
Q.VV32_2d = dim_scale * (sum((E / kT0) .* RVV32_2d) + RVV32_2d_partner_E / kT0);
Q.VV34_2d = dim_scale * (sum((E / kT0) .* RVV34_2d) + RVV34_2d_partner_E / kT0);
Q.VV32_4d = dim_scale * (sum((E / kT0) .* RVV32_4d) + RVV32_4d_partner_E / kT0);
end

function kb = k_vt_bwd(kf, i, f, T, AD)
kb = kf * (AD.stw(i) / AD.stw(f)) * exp((AD.e1234(f) - AD.e1234(i)) / (AD.k * T));
end

function test_relaxation_terms(varargin)
%TEST_RELAXATION_TERMS Print diagnostic relaxation terms for CH4 STS models.
%
% Run the full diagnostic:
%   test_relaxation_terms
%
% Run a fast smoke check on a truncated state space:
%   test_relaxation_terms('quick')
%
% The source assembly below mirrors the VT2, VT4, and VV34d blocks in
% rpart_mt_sts.m and rpart_mt_3t_sts.m. It intentionally prints only the
% general state-specific moments; reduced/regression formulas are not used.

quick_mode = any(strcmpi(varargin, 'quick'));

this_dir = fileparts(mfilename('fullpath'));
addpath(this_dir);
addpath(fullfile(this_dir, 'fho_model'));

AD = build_test_data();

if quick_mode
    AD = truncate_state_space(AD, [2, 3, 2, 4]);
    n_steps_vt = 400;
    n_steps_vv = 300;
    fprintf('Running quick smoke diagnostic on truncated state space.\n');
else
    n_steps_vt = 8000;
    n_steps_vv = 6000;
    fprintf('Running full diagnostic on input_data state space.\n');
end

fprintf('T0 = %.6g K, p0 = %.6g Pa, n0 = %.6e m^-3, tau = %.6e s\n', ...
    AD.T0, AD.p0, AD.n0, AD.tau);
fprintf('FHO params: alpha = %.6e m^-1, e_m = %.6e K, steric2 = %.6g, steric4 = %.6g, steric_vv_34d = %.6g\n', ...
    AD.fho_alpha, AD.fho_e_m, AD.fho_steric2, AD.fho_steric4, AD.fho_steric_vv_34d);
fprintf('Rate quadrature steps: VT = %d, VV = %d\n\n', n_steps_vt, n_steps_vv);

n_report = 1.0e20;
tau_report = mean_collision_time(n_report, AD.T0, AD);
dimensional_factor = n_report * AD.k * AD.T0 / tau_report;
nondim_factor = 1 / dimensional_factor;
fprintf('Printed values are nondimensional RHS moments used by the MATLAB MT equations:\n');
fprintf('  x_i = normalized multi-temperature Boltzmann population.\n');
fprintf('  r_i^nd = (n0*tau) * sum_pair_fluxes(x,k), where k is in m^3/s.\n');
fprintf('  R_13^nd = sum_i ((i1*eps1+i3*eps3)/(k*T0)) * r_i^nd.\n');
fprintf('  R_24^nd = sum_i ((i2*eps2+i4*eps4)/(k*T0)) * r_i^nd.\n');
fprintf('  These are the quantities printed below.\n\n');
fprintf('If Python gives dimensional energy source Q_dim [J/(m^3*s)] at n_ref, use:\n');
fprintf('  Q_nd = Q_dim * tau_ref/(n_ref*k*T0), tau_ref = 1/(4*n_ref*pi*r0^2*sqrt(k*T0/(pi*m))).\n');
fprintf('  For n_ref = %.3e m^-3: tau_ref = %.6e s, Q_nd = Q_dim * %.6e.\n', ...
    n_report, tau_report, nondim_factor);
fprintf('  Equivalently, Q_dim = Q_nd * %.6e J/(m^3*s).\n\n', dimensional_factor);

test_cases_3t = [
    800,  700,  600;
    500,  800,  400;
    1200, 800, 1400;
];

test_cases_2t = [
    800,  600;
    500,  800;
    1200, 1000;
];

eps = mode_energies(AD);
expected_ratio = -2 * eps.eps4 / eps.eps3;

fprintf('%s\n', repmat('=', 1, 128));
fprintf('3T MODEL: selected general nondimensional moments\n');
fprintf('%s\n', repmat('=', 1, 128));
fprintf('%6s %6s %6s %14s %10s %14s %10s %14s %14s %10s\n', ...
    'T', 'T13', 'T24', 'VT2 R24nd', 'VT2 R13', 'VT4 R24nd', 'VT4 R13', ...
    'VV34d R13nd', 'VV34d R24nd', 'ratio err');
fprintf('%s\n', repmat('-', 1, 128));

for icase = 1:size(test_cases_3t, 1)
    T = test_cases_3t(icase, 1);
    T13 = test_cases_3t(icase, 2);
    T24 = test_cases_3t(icase, 3);

    terms = compute_terms_3t(T, T13, T24, AD, n_steps_vt, n_steps_vv);
    ratio_err = relative_error(safe_ratio(terms.VV34d.R24, terms.VV34d.R13), expected_ratio);

    fprintf('%6.0f %6.0f %6.0f %14.6e %10.2e %14.6e %10.2e %14.6e %14.6e %10.2e\n', ...
        T, T13, T24, ...
        terms.VT2.R24, terms.VT2.R13, ...
        terms.VT4.R24, terms.VT4.R13, ...
        terms.VV34d.R13, terms.VV34d.R24, ratio_err);

    assert_selected_terms(terms, expected_ratio, sprintf('3T case %d', icase));
end

fprintf('\n%s\n', repmat('=', 1, 128));
fprintf('2T MODEL: selected general nondimensional moments\n');
fprintf('%s\n', repmat('=', 1, 128));
fprintf('%6s %6s %14s %10s %14s %10s %14s %14s %14s %10s\n', ...
    'T', 'Tv', 'VT2 Rvibnd', 'VT2 off', 'VT4 Rvibnd', 'VT4 off', ...
    'VV34d Rvibnd', 'VV34d R3nd', 'VV34d R4nd', 'ratio err');
fprintf('%s\n', repmat('-', 1, 128));

for icase = 1:size(test_cases_2t, 1)
    T = test_cases_2t(icase, 1);
    Tv = test_cases_2t(icase, 2);

    terms = compute_terms_2t(T, Tv, AD, n_steps_vt, n_steps_vv);
    vt2_off = max(abs([terms.VT2.R1, terms.VT2.R3, terms.VT2.R4]));
    vt4_off = max(abs([terms.VT4.R1, terms.VT4.R2, terms.VT4.R3]));
    ratio_err = relative_error(safe_ratio(terms.VV34d.R4, terms.VV34d.R3), expected_ratio);

    fprintf('%6.0f %6.0f %14.6e %10.2e %14.6e %10.2e %14.6e %14.6e %14.6e %10.2e\n', ...
        T, Tv, ...
        terms.VT2.Rvibr, vt2_off, ...
        terms.VT4.Rvibr, vt4_off, ...
        terms.VV34d.Rvibr, terms.VV34d.R3, terms.VV34d.R4, ratio_err);

    assert_selected_terms(terms, expected_ratio, sprintf('2T case %d', icase));
end

fprintf('%s\n', repmat('=', 1, 128));
fprintf('Expected VV34d ratio R4/R3 = R24/R13 = %.12e\n', expected_ratio);
fprintf('All selected relaxation-term checks passed.\n');

end

function AD = build_test_data()
AD = input_data();

AD.T0 = 1000;
AD.p0 = 101325;
AD.n0 = AD.p0 / (AD.k * AD.T0);
AD.tau = mean_collision_time(AD.n0, AD.T0, AD);

AD.sw_rt = 'fho';
AD.fho_alpha = 5.174e10;
AD.fho_e_m = 776.4;
AD.fho_steric2 = 0.002612;
AD.fho_steric4 = 0.028546;
AD.fho_steric_vv_34s = 1.0;
AD.fho_steric_vv_34d = 0.24;
AD.fho_gamma_vv34d = 1.0;
AD.fho_steric_vv_34_4d = 1.0;

AD = ch4_relax_topology(AD);
end

function tau = mean_collision_time(n, T, AD)
sigma0 = pi * AD.r0^2;
tau = (4 * n * sigma0 * sqrt(AD.k * T / (pi * AD.m)))^(-1);
end

function AD = truncate_state_space(AD, limits)
mask = AD.inds(:, 1) < limits(1) & ...
       AD.inds(:, 2) < limits(2) & ...
       AD.inds(:, 3) < limits(3) & ...
       AD.inds(:, 4) < limits(4);

AD.inds = AD.inds(mask, :);
AD.e1234 = AD.e1234(mask);
AD.stw = AD.stw(mask);
AD.lch4 = limits;
AD.klop = nnz(mask);

remove_fields = {'indj_vt2', 'indl_vt4', 'indvv34', 'indvv34d'};
for ifield = 1:numel(remove_fields)
    if isfield(AD, remove_fields{ifield})
        AD = rmfield(AD, remove_fields{ifield});
    end
end
AD = ch4_relax_topology(AD);
end

function terms = compute_terms_2t(T, Tv, AD, n_steps_vt, n_steps_vv)
x = boltzmann_2t(Tv, AD);
sources = selected_sources(T, x, AD, n_steps_vt, n_steps_vv);
terms = summarize_terms(sources, AD);
terms.model = '2T';
terms.T = T;
terms.Tv = Tv;
end

function terms = compute_terms_3t(T, T13, T24, AD, n_steps_vt, n_steps_vv)
x = boltzmann_3t(T13, T24, AD);
sources = selected_sources(T, x, AD, n_steps_vt, n_steps_vv);
terms = summarize_terms(sources, AD);
terms.model = '3T';
terms.T = T;
terms.T13 = T13;
terms.T24 = T24;
end

function x = boltzmann_2t(Tv, AD)
E = AD.e1234(:);
f = AD.stw(:) .* exp(-E ./ (AD.k * Tv));
x = f ./ sum(f);
end

function x = boltzmann_3t(T13, T24, AD)
eps = mode_energies(AD);
I0 = AD.inds(:, 1);
J0 = AD.inds(:, 2);
K0 = AD.inds(:, 3);
L0 = AD.inds(:, 4);

E13 = I0 * eps.eps1 + K0 * eps.eps3;
E24 = J0 * eps.eps2 + L0 * eps.eps4;

f = AD.stw(:) .* exp(-(E13 ./ (AD.k * T13) + E24 ./ (AD.k * T24)));
x = f ./ sum(f);
end

function sources = selected_sources(T, x, AD, n_steps_vt, n_steps_vv)
n_scale = AD.n0 * AD.tau;

sources.VT2 = vt_source(T, x, AD, AD.indj_vt2, 2, AD.fho_steric2, n_scale, n_steps_vt);
sources.VT4 = vt_source(T, x, AD, AD.indl_vt4, 4, AD.fho_steric4, n_scale, n_steps_vt);
sources.VV34d = vv34d_source(T, x, AD, n_scale, n_steps_vv);
end

function R = vt_source(T, x, AD, index_next, mode, steric, n_scale, n_steps)
N = numel(AD.e1234);
R = zeros(N, 1);
src = find(~isnan(index_next));

for ii = 1:numel(src)
    r = src(ii);
    dst = index_next(r);
    q_high = AD.inds(dst, mode);

    kf = cached_rate_vt(T, q_high, mode, steric, AD.fho_alpha, AD.fho_e_m, n_steps) * n_scale;
    kr = detailed_balance_state(kf, dst, r, T, AD);

    R(r) = R(r) + x(dst) * kf - x(r) * kr;
    R(dst) = R(dst) + x(r) * kr - x(dst) * kf;
end
end

function R = vv34d_source(T, x, AD, n_scale, n_steps)
N = numel(AD.e1234);
R = zeros(N, 1);
src = find(~isnan(AD.indvv34d));
sk0 = [0, 0, 0, 0];

for ii = 1:numel(src)
    r = src(ii);
    dst = AD.indvv34d(r);
    si = AD.inds(r, :);
    sf = AD.inds(dst, :);

    kf = cached_rate_vv34d(T, si, sf, sk0, AD.fho_steric_vv_34d, ...
        AD.fho_alpha, AD.fho_e_m, AD.fho_gamma_vv34d, n_steps) * n_scale;
    kr = detailed_balance_state(kf, r, dst, T, AD);

    R(r) = R(r) + x(dst) * kr - x(r) * kf;
    R(dst) = R(dst) + x(r) * kf - x(dst) * kr;
end
end

function terms = summarize_terms(sources, AD)
[W1, W2, W3, W4, Wvibr] = energy_weights(AD);

terms.VT2 = summarize_one_source(sources.VT2, W1, W2, W3, W4, Wvibr);
terms.VT4 = summarize_one_source(sources.VT4, W1, W2, W3, W4, Wvibr);
terms.VV34d = summarize_one_source(sources.VV34d, W1, W2, W3, W4, Wvibr);
end

function one = summarize_one_source(R, W1, W2, W3, W4, Wvibr)
one.R1 = sum(W1 .* R);
one.R2 = sum(W2 .* R);
one.R3 = sum(W3 .* R);
one.R4 = sum(W4 .* R);
one.R13 = one.R1 + one.R3;
one.R24 = one.R2 + one.R4;
one.Rvibr = sum(Wvibr .* R);
one.Rvibr_modes = one.R1 + one.R2 + one.R3 + one.R4;
one.population_balance = sum(R);
end

function [W1, W2, W3, W4, Wvibr] = energy_weights(AD)
eps = mode_energies(AD);
kT0 = AD.k * AD.T0;

I0 = AD.inds(:, 1);
J0 = AD.inds(:, 2);
K0 = AD.inds(:, 3);
L0 = AD.inds(:, 4);

W1 = I0 * eps.eps1 / kT0;
W2 = J0 * eps.eps2 / kT0;
W3 = K0 * eps.eps3 / kT0;
W4 = L0 * eps.eps4 / kT0;
Wvibr = AD.e1234(:) / kT0;
end

function eps = mode_energies(AD)
eps.eps1 = AD.e1000 - AD.e0000;
eps.eps2 = AD.e0100 - AD.e0000;
eps.eps3 = AD.e0010 - AD.e0000;
eps.eps4 = AD.e0001 - AD.e0000;
end

function k = cached_rate_vt(T, q_high, mode, steric, alpha, e_m, n_steps)
persistent cache
if isempty(cache)
    cache = containers.Map('KeyType', 'char', 'ValueType', 'double');
end

key = sprintf('VT_%d_%.12g_%d_%.12g_%.12g_%.12g_%d', ...
    mode, T, q_high, steric, alpha, e_m, n_steps);

if isKey(cache, key)
    k = cache(key);
else
    k = rates_fho_vt(T, q_high, q_high - 1, mode, steric, alpha, e_m, ...
        'n_steps', n_steps);
    cache(key) = k;
end
end

function k = cached_rate_vv34d(T, si, sf, sk0, steric, alpha, e_m, gamma, n_steps)
persistent cache
if isempty(cache)
    cache = containers.Map('KeyType', 'char', 'ValueType', 'double');
end

key = sprintf('VV34d_%.12g_%d_%d_%.12g_%.12g_%.12g_%.12g_%d', ...
    T, si(3), si(4), steric, alpha, e_m, gamma, n_steps);

if isKey(cache, key)
    k = cache(key);
else
    k = rates_fho_vv(T, si, sf, sk0, sk0, steric, alpha, e_m, ...
        'n_steps', n_steps, 'gamma', gamma);
    cache(key) = k;
end
end

function kb = detailed_balance_state(kf, i, f, T, AD)
kb = kf * (AD.stw(i) / AD.stw(f)) * exp((AD.e1234(f) - AD.e1234(i)) / (AD.k * T));
end

function assert_selected_terms(terms, expected_ratio, label)
zero_tol = 1e-9;
ratio_tol = 1e-8;
balance_tol = 1e-10;
tiny = realmin;

assert(isfinite_source(terms.VT2), '%s VT2 contains a non-finite value.', label);
assert(isfinite_source(terms.VT4), '%s VT4 contains a non-finite value.', label);
assert(isfinite_source(terms.VV34d), '%s VV34d contains a non-finite value.', label);

assert(abs(terms.VT2.R13) <= zero_tol * max(abs(terms.VT2.R24), 1), ...
    '%s VT2 should not contribute to R13.', label);
assert(abs(terms.VT4.R13) <= zero_tol * max(abs(terms.VT4.R24), 1), ...
    '%s VT4 should not contribute to R13.', label);
assert(max(abs([terms.VV34d.R1, terms.VV34d.R2])) <= zero_tol * max(abs(terms.VV34d.R3) + abs(terms.VV34d.R4), 1), ...
    '%s VV34d should not contribute to modes 1 or 2.', label);

assert(abs(terms.VT2.population_balance) <= balance_tol * max(norm_source(terms.VT2), 1), ...
    '%s VT2 does not conserve population.', label);
assert(abs(terms.VT4.population_balance) <= balance_tol * max(norm_source(terms.VT4), 1), ...
    '%s VT4 does not conserve population.', label);
assert(abs(terms.VV34d.population_balance) <= balance_tol * max(norm_source(terms.VV34d), 1), ...
    '%s VV34d does not conserve population.', label);

if abs(terms.VV34d.R13) > tiny
    ratio_3t = terms.VV34d.R24 / terms.VV34d.R13;
    assert(relative_error(ratio_3t, expected_ratio) <= ratio_tol, ...
        '%s VV34d R24/R13 ratio is inconsistent.', label);
end
if abs(terms.VV34d.R3) > tiny
    ratio_2t = terms.VV34d.R4 / terms.VV34d.R3;
    assert(relative_error(ratio_2t, expected_ratio) <= ratio_tol, ...
        '%s VV34d R4/R3 ratio is inconsistent.', label);
end

assert(abs(terms.VT2.Rvibr - terms.VT2.Rvibr_modes) <= zero_tol * max(abs(terms.VT2.Rvibr), 1), ...
    '%s VT2 modal sum does not match vibrational moment.', label);
assert(abs(terms.VT4.Rvibr - terms.VT4.Rvibr_modes) <= zero_tol * max(abs(terms.VT4.Rvibr), 1), ...
    '%s VT4 modal sum does not match vibrational moment.', label);
assert(abs(terms.VV34d.Rvibr - terms.VV34d.Rvibr_modes) <= zero_tol * max(abs(terms.VV34d.Rvibr), 1), ...
    '%s VV34d modal sum does not match vibrational moment.', label);
end

function tf = isfinite_source(one)
vals = [one.R1, one.R2, one.R3, one.R4, one.R13, one.R24, ...
    one.Rvibr, one.Rvibr_modes, one.population_balance];
tf = all(isfinite(vals));
end

function nrm = norm_source(one)
vals = [one.R1, one.R2, one.R3, one.R4, one.R13, one.R24, one.Rvibr];
nrm = norm(vals, inf);
end

function r = safe_ratio(num, den)
if abs(den) <= realmin
    r = NaN;
else
    r = num / den;
end
end

function err = relative_error(value, reference)
if ~isfinite(value)
    err = NaN;
else
    err = abs(value - reference) / max(abs(reference), realmin);
end
end

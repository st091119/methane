function test_ch4_relax_terms_simplified(varargin)
%TEST_CH4_RELAX_TERMS_SIMPLIFIED Print reduced 3T relaxation-term diagnostics.
%
% Run:
%   test_ch4_relax_terms_simplified
%
% Faster smoke run:
%   test_ch4_relax_terms_simplified('quick')
%
% The values printed here are the simplified VT2, VT4, VV34s, VV34d, and
% VV34_4d terms from ch4_relax_terms_simplified.m. They are intended to be
% compared directly with the analogous notebook functions after applying
% the same nondimensionalization:
%
%   R_dim = n0^2 * eps_mode * reduced_sum
%   R_nd  = R_dim * tau/(n0*k*T0)
%         = n0*tau*eps_mode*reduced_sum/(k*T0).

quick_mode = any(strcmpi(varargin, 'quick'));

this_dir = fileparts(mfilename('fullpath'));
addpath(this_dir);
addpath(fullfile(this_dir, 'fho_model'));

AD = build_test_data();
if quick_mode
    AD = truncate_state_space(AD, [2, 2, 3, 4]);
    n_steps_vv = 300;
    fprintf('Running quick simplified-term smoke diagnostic.\n');
else
    n_steps_vv = 6000;
    fprintf('Running simplified-term diagnostic on input_data state space.\n');
end

ch4_relax_terms_simplified('clear_cache');

eps = mode_energies_local(AD);
scale_nd = AD.n0 * AD.tau / (AD.k * AD.T0);
scale_dim = AD.n0^2;

fprintf('T0 = %.6g K, p0 = %.6g Pa, n0 = %.6e m^-3, tau = %.6e s\n', ...
    AD.T0, AD.p0, AD.n0, AD.tau);
fprintf('Nondimensionalization: R_nd = R_dim*tau/(n0*k*T0) = %.6e * eps * reduced_sum\n', ...
    scale_nd);
fprintf('Dimensionalization:    R_dim = n0^2 * eps * reduced_sum = %.6e * eps * reduced_sum\n', ...
    scale_dim);
fprintf(['FHO VV params: alpha = %.6e m^-1, e_m = %.6e K, ' ...
    'steric34s = %.6g, steric34d = %.6g, gamma34d = %.6g, steric34_4d = %.6g\n'], ...
    AD.fho_alpha, AD.fho_e_m, AD.fho_steric_vv_34s, AD.fho_steric_vv_34d, ...
    AD.fho_gamma_vv34d, AD.fho_steric_vv_34_4d);
fprintf('VV quadrature steps = %d\n\n', n_steps_vv);

test_cases_3t = [
    800,  700,  600;
    500,  800,  400;
    1200, 800, 1400;
];

expected_ratio_s = -eps.eps4 / eps.eps3;
expected_ratio_d = -2 * eps.eps4 / eps.eps3;

fprintf('%s\n', repmat('=', 1, 128));
fprintf('3T SIMPLIFIED RELAXATION TERMS: reduced sums and nondimensional RHS moments\n');
fprintf('%s\n', repmat('=', 1, 128));
fprintf('%6s %6s %6s %13s %13s %13s %13s %13s %13s %13s %13s\n', ...
    'T', 'T13', 'T24', 'S2', 'VT2 R24', 'S4', 'VT4 R24', ...
    'A34s', 'VV34s R13', 'VV34s R24', 'ratio err');
fprintf('%s\n', repmat('-', 1, 128));

for icase = 1:size(test_cases_3t, 1)
    T = test_cases_3t(icase, 1);
    T13 = test_cases_3t(icase, 2);
    T24 = test_cases_3t(icase, 3);

    terms = ch4_relax_terms_simplified(T, T13, T24, AD, 'n_steps_vv', n_steps_vv);
    ratio_err_s = relative_error(safe_ratio(terms.VV34s.R24, terms.VV34s.R13), expected_ratio_s);

    fprintf('%6.0f %6.0f %6.0f %13.6e %13.6e %13.6e %13.6e %13.6e %13.6e %13.6e %13.2e\n', ...
        T, T13, T24, ...
        terms.VT2.reduced_sum, terms.VT2.R24, ...
        terms.VT4.reduced_sum, terms.VT4.R24, ...
        terms.VV34s.reduced_sum, terms.VV34s.R13, terms.VV34s.R24, ratio_err_s);

    assert_terms(terms, expected_ratio_s, expected_ratio_d, sprintf('case %d', icase));
end

fprintf('\n%s\n', repmat('=', 1, 128));
fprintf('3T SIMPLIFIED VV DEFECT TERMS\n');
fprintf('%s\n', repmat('=', 1, 128));
fprintf('%6s %6s %6s %13s %13s %13s %13s %13s %13s %13s\n', ...
    'T', 'T13', 'T24', 'A34d', 'VV34d R13', 'VV34d R24', ...
    'B34_4d', 'VV34_4d R13', 'VV34_4d R24', 'ratio err');
fprintf('%s\n', repmat('-', 1, 128));

for icase = 1:size(test_cases_3t, 1)
    T = test_cases_3t(icase, 1);
    T13 = test_cases_3t(icase, 2);
    T24 = test_cases_3t(icase, 3);

    terms = ch4_relax_terms_simplified(T, T13, T24, AD, 'n_steps_vv', n_steps_vv);
    ratio_err_d = relative_error(safe_ratio(terms.VV34d.R24, terms.VV34d.R13), expected_ratio_d);
    ratio_err_b = relative_error(safe_ratio(terms.VV34_4d.R24, terms.VV34_4d.R13), expected_ratio_d);

    fprintf('%6.0f %6.0f %6.0f %13.6e %13.6e %13.6e %13.6e %13.6e %13.6e %13.2e\n', ...
        T, T13, T24, ...
        terms.VV34d.reduced_sum, terms.VV34d.R13, terms.VV34d.R24, ...
        terms.VV34_4d.reduced_sum, terms.VV34_4d.R13, terms.VV34_4d.R24, ...
        max(ratio_err_d, ratio_err_b));
end

fprintf('\n%s\n', repmat('=', 1, 128));
fprintf('3T TOTALS FROM SIMPLIFIED TERMS\n');
fprintf('%s\n', repmat('=', 1, 128));
fprintf('%6s %6s %6s %14s %14s %14s %14s\n', ...
    'T', 'T13', 'T24', 'R13 nd', 'R24 nd', 'R13 dim', 'R24 dim');
fprintf('%s\n', repmat('-', 1, 128));

for icase = 1:size(test_cases_3t, 1)
    T = test_cases_3t(icase, 1);
    T13 = test_cases_3t(icase, 2);
    T24 = test_cases_3t(icase, 3);

    terms = ch4_relax_terms_simplified(T, T13, T24, AD, 'n_steps_vv', n_steps_vv);
    fprintf('%6.0f %6.0f %6.0f %14.6e %14.6e %14.6e %14.6e\n', ...
        T, T13, T24, terms.R13, terms.R24, terms.R13_dim, terms.R24_dim);
end

fprintf('%s\n', repmat('=', 1, 128));
fprintf('Expected VV34s ratio R24/R13 = %.12e\n', expected_ratio_s);
fprintf('Expected VV34d and VV34_4d ratio R24/R13 = %.12e\n', expected_ratio_d);
fprintf('All simplified 3T relaxation-term checks passed.\n');

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
AD.fho_steric_vv_32_2d = 1.0;
AD.fho_steric_vv_34_2d = 1.0;
AD.fho_steric_vv_32_4d = 1.0;
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
end

function eps = mode_energies_local(AD)
eps.eps1 = AD.e1000 - AD.e0000;
eps.eps2 = AD.e0100 - AD.e0000;
eps.eps3 = AD.e0010 - AD.e0000;
eps.eps4 = AD.e0001 - AD.e0000;
end

function assert_terms(terms, expected_ratio_s, expected_ratio_d, label)
assert(isfinite(terms.VT2.R24), '%s VT2 R24 is non-finite.', label);
assert(isfinite(terms.VT4.R24), '%s VT4 R24 is non-finite.', label);
assert(isfinite(terms.VV34s.R13) && isfinite(terms.VV34s.R24), '%s VV34s is non-finite.', label);
assert(isfinite(terms.VV34d.R13) && isfinite(terms.VV34d.R24), '%s VV34d is non-finite.', label);
assert(isfinite(terms.VV34_4d.R13) && isfinite(terms.VV34_4d.R24), '%s VV34_4d is non-finite.', label);
assert(abs(terms.VT2.R13) <= 10 * eps, '%s VT2 should not contribute to R13.', label);
assert(abs(terms.VT4.R13) <= 10 * eps, '%s VT4 should not contribute to R13.', label);

check_ratio(terms.VV34s.R24, terms.VV34s.R13, expected_ratio_s, [label ' VV34s']);
check_ratio(terms.VV34d.R24, terms.VV34d.R13, expected_ratio_d, [label ' VV34d']);
check_ratio(terms.VV34_4d.R24, terms.VV34_4d.R13, expected_ratio_d, [label ' VV34_4d']);
end

function check_ratio(num, den, expected, label)
if abs(den) <= 1e-300
    return
end
actual = num / den;
assert(relative_error(actual, expected) <= 1e-11, '%s ratio is inconsistent.', label);
end

function value = safe_ratio(num, den)
if abs(den) <= 1e-300
    value = NaN;
else
    value = num / den;
end
end

function err = relative_error(actual, expected)
if ~isfinite(actual)
    err = NaN;
else
    err = abs(actual - expected) / max(abs(expected), 1e-300);
end
end

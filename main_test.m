
function results = main_test(varargin)
%MAIN_TEST Fast cumulative check of simplified CH4 2T/3T relaxation models.
%
% Run the fast diagnostic:
%   results = main_test;
%
% Run a heavier diagnostic:
%   results = main_test('full');
%
% Useful name-value options:
%   't_fin'       final physical time [s]
%   'n_steps_vt'  FHO quadrature steps for direct VT averaged sums
%   'n_steps_vv'  FHO quadrature steps for VV reduced sums
%   'make_plots'  true/false
%   'n_output'    number of log-spaced output times
%   'plot_t_min'  first positive time shown on log plots [s]
%
% The simplified runs use ch4_relax_terms_simplified.m with a selected
% process list, so VT-only cases do not evaluate the VV sums.

cfg = parse_options(varargin{:});

this_dir = fileparts(mfilename('fullpath'));
addpath(this_dir);
addpath(fullfile(this_dir, 'fho_model'));

AD = build_test_data(cfg);
ch4_relax_terms_simplified('clear_cache');

tspan = output_time_grid(cfg) / AD.tau;
Y0_2t = [1; cfg.Tv0 / cfg.T0];
Y0_3t = [1; cfg.Tv0 / cfg.T0; cfg.Tv0 / cfg.T0];

ode_options = odeset( ...
    'RelTol', cfg.rel_tol, ...
    'AbsTol', cfg.abs_tol, ...
    'MaxStep', cfg.max_step / AD.tau);

run_cases = make_run_cases();
results = repmat(empty_result(), size(run_cases));
T_eq = equilibrium_temperature(AD, cfg.T0, cfg.Tv0);

fprintf('main_test: %s mode\n', cfg.mode_name);
fprintf('Spectrum: %s, lch4 = [%s], klop = %d\n', AD.spectrum_mode, num2str(AD.lch4), AD.klop);
fprintf('T0 = %.6g K, Tv0 = %.6g K, p0 = %.6g Pa, t_fin = %.6e s\n', ...
    cfg.T0, cfg.Tv0, cfg.p0, cfg.t_fin);
fprintf('Equilibrium temperature from total energy conservation: %.6f K\n', T_eq);
fprintf('n0 = %.6e m^-3, tau = %.6e s\n', AD.n0, AD.tau);
fprintf('Rate quadrature steps: VT = %d, VV = %d\n\n', cfg.n_steps_vt, cfg.n_steps_vv);
print_vt_source_diagnostic(AD, cfg);

for icase = 1:numel(run_cases)
    c = run_cases(icase);
    fprintf('--- Solving %s ---\n', c.name);
    t_start = tic;

    switch c.kind
        case 'lt_wang'
            AD_run = AD;
            AD_run.sw_rt = 'vt_rel_time_wang';
            [X, Y] = ode15s(@(t, y) rpart_mt_lt(t, y, AD_run), tspan, Y0_2t, ode_options);
        case 'lt_fho'
            AD_run = AD;
            AD_run.sw_rt = 'fho';
            [X, Y] = ode15s(@(t, y) rpart_mt_lt(t, y, AD_run), tspan, Y0_2t, ode_options);
        case 'simplified_2t'
            [X, Y] = ode15s(@(t, y) rhs_simplified_2t(t, y, AD, c.processes, cfg), ...
                tspan, Y0_2t, ode_options);
        case 'simplified_3t'
            [X, Y] = ode15s(@(t, y) rhs_simplified_3t(t, y, AD, c.processes, cfg), ...
                tspan, Y0_3t, ode_options);
        otherwise
            error('main_test:unknown_case_kind', 'Unknown case kind: %s', c.kind);
    end

    results(icase) = pack_result(c, X, Y, AD, toc(t_start));
    print_result_summary(results(icase), T_eq);
end

print_equilibrium_table(results, T_eq);

if cfg.make_plots
    plot_main_style(results, T_eq, cfg);
    plot_split_style(results, T_eq, cfg, 'All cumulative cases');
end

fprintf('\nmain_test finished.\n');

end

function cfg = parse_options(varargin)
cfg = struct();
cfg.mode_name = 'fast';
cfg.fast = true;
cfg.T0 = 800;
cfg.Tv0 = 300;
cfg.p0 = 101325;
cfg.t_fin = 1.0e-5;
cfg.n_steps_vt = 300;
cfg.n_steps_vv = 100;
cfg.rel_tol = 1e-6;
cfg.abs_tol = 1e-8;
cfg.max_step = 2.0e-4;
cfg.n_output = 360;
cfg.plot_t_min = 1.0e-16;
cfg.make_plots = true;
cfg.truncate_limits = [4, 4, 4, 4];

i = 1;
while i <= numel(varargin)
    arg = varargin{i};
    if ischar(arg) || isstring(arg)
        name = lower(char(string(arg)));
        switch name
            case {'fast', 'quick'}
                cfg.mode_name = 'fast';
                cfg.fast = true;
                i = i + 1;
                continue
            case 'full'
                cfg.mode_name = 'full';
                cfg.fast = false;
                cfg.n_steps_vt = 8000;
                cfg.n_steps_vv = 6000;
                cfg.rel_tol = 1e-8;
                cfg.abs_tol = 1e-8;
                cfg.max_step = 2.0e-4;
                i = i + 1;
                continue
            case 't_fin'
                cfg.t_fin = varargin{i + 1};
            case 'n_steps_vv'
                cfg.n_steps_vv = varargin{i + 1};
            case {'n_steps_vt', 'n_steps_vt_check'}
                cfg.n_steps_vt = varargin{i + 1};
            case 'make_plots'
                cfg.make_plots = logical(varargin{i + 1});
            case 'n_output'
                cfg.n_output = varargin{i + 1};
            case 'plot_t_min'
                cfg.plot_t_min = varargin{i + 1};
            case 't0'
                cfg.T0 = varargin{i + 1};
            case 'tv0'
                cfg.Tv0 = varargin{i + 1};
            case 'p0'
                cfg.p0 = varargin{i + 1};
            case 'rel_tol'
                cfg.rel_tol = varargin{i + 1};
            case 'abs_tol'
                cfg.abs_tol = varargin{i + 1};
            case 'max_step'
                cfg.max_step = varargin{i + 1};
            case 'truncate_limits'
                cfg.truncate_limits = varargin{i + 1};
            otherwise
                error('main_test:unknown_option', 'Unknown option: %s', name);
        end
        i = i + 2;
    else
        error('main_test:bad_option', 'Options must be strings or name-value pairs.');
    end
end
end

function t = output_time_grid(cfg)
t_min = max(cfg.plot_t_min, realmin);
if t_min >= cfg.t_fin
    t = [0, cfg.t_fin];
    return
end
t = [0, logspace(log10(t_min), log10(cfg.t_fin), cfg.n_output)];
t = unique(t, 'stable');
end

function AD = build_test_data(cfg)
AD = input_data();

if cfg.fast
    limits = min(AD.lch4, cfg.truncate_limits);
    if any(AD.lch4 ~= limits)
        AD = truncate_state_space(AD, limits);
        AD.spectrum_mode = [AD.spectrum_mode '_truncated_for_main_test'];
    end
end

AD.T0 = cfg.T0;
AD.p0 = cfg.p0;
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

function T_eq = equilibrium_temperature(AD, T0, Tv0)
target = 3 * AD.k * T0 + mean_vibrational_energy(Tv0, AD);
energy_residual = @(T) 3 * AD.k * T + mean_vibrational_energy(T, AD) - target;

lo = 1;
hi = max(T0, Tv0) * 2;
while energy_residual(hi) < 0
    hi = hi * 2;
end
T_eq = fzero(energy_residual, [lo, hi]);
end

function Ev = mean_vibrational_energy(T, AD)
w = AD.stw(:);
E = AD.e1234(:);
xi = -E ./ (AD.k * T);
x = normalized_boltzmann(w, xi);
Ev = sum(x .* E);
end

function run_cases = make_run_cases()
run_cases = struct( ...
    'name', { ...
        'LT-FHO 2T reference', ...
        '2T STS-averaged VT2+VT4 direct', ...
        '3T VT2+VT4+VV34s', ...
        '3T VT2+VT4+VV34s+VV34d', ...
        '3T VT2+VT4+VV34s+VV34d+VV34_4d'}, ...
    'kind', { ...
        'lt_fho', ...
        'simplified_2t', ...
        'simplified_3t', ...
        'simplified_3t', ...
        'simplified_3t'}, ...
    'dim', {2, 2, 3, 3, 3}, ...
    'processes', { ...
        {{}}, ...
        {{'VT2', 'VT4'}}, ...
        {{'VT2', 'VT4', 'VV34s'}}, ...
        {{'VT2', 'VT4', 'VV34s', 'VV34d'}}, ...
        {{'VT2', 'VT4', 'VV34s', 'VV34d', 'VV34_4d'}}} ...
);
end

function print_vt_source_diagnostic(AD, cfg)
T = cfg.T0;
Tv = cfg.Tv0;

reg = ch4_relax_terms_simplified(T, Tv, AD, 'processes', {'VT2', 'VT4'});
[S2_direct, S4_direct] = direct_vt_sums(T, Tv, AD, cfg.n_steps_vt);

eps = mode_energies(AD);
scale_nd = AD.n0 * AD.tau / (AD.k * AD.T0);
R2_direct = scale_nd * eps.eps2 * S2_direct;
R4_direct = scale_nd * eps.eps4 * S4_direct;
Rv_direct = R2_direct + R4_direct;

fprintf('VT source check at T = %.1f K, Tv = %.1f K:\n', T, Tv);
fprintf('  S2 regression/direct = %.6e / %.6e, ratio = %.3e\n', ...
    reg.VT2.reduced_sum, S2_direct, safe_ratio(reg.VT2.reduced_sum, S2_direct));
fprintf('  S4 regression/direct = %.6e / %.6e, ratio = %.3e\n', ...
    reg.VT4.reduced_sum, S4_direct, safe_ratio(reg.VT4.reduced_sum, S4_direct));
fprintf('  Rvibr_nd regression/direct = %.6e / %.6e, ratio = %.3e\n\n', ...
    reg.Rvibr, Rv_direct, safe_ratio(reg.Rvibr, Rv_direct));
fprintf('  main_test uses the direct STS-averaged VT source in the ODE.\n\n');
end

function [S2, S4] = direct_vt_sums(T, Tstar, AD, n_steps)
eps = mode_energies(AD);
S2 = direct_vt_sum_mode(T, Tstar, 2, eps.eps2, AD, n_steps);
S4 = direct_vt_sum_mode(T, Tstar, 4, eps.eps4, AD, n_steps);
end

function S = direct_vt_sum_mode(T, Tstar, mode, eps_m, AD, n_steps)
switch mode
    case 2
        level_max = AD.lch4(2) - 1;
        steric = AD.fho_steric2;
    case 4
        level_max = AD.lch4(4) - 1;
        steric = AD.fho_steric4;
    otherwise
        error('main_test:bad_vt_mode', 'mode must be 2 or 4.');
end

Z = 0;
acc = 0;
for q = 0:level_max
    s_q = stat_weight_mode(q, mode);
    boltz = exp(-q * eps_m / (AD.k * Tstar));
    k_down = 0;
    if q > 0
        k_down = rates_fho_vt(T, q, q - 1, mode, steric, AD.fho_alpha, AD.fho_e_m, ...
            'n_steps', n_steps);
    end
    k_up = 0;
    if q < level_max
        k_down_next = rates_fho_vt(T, q + 1, q, mode, steric, AD.fho_alpha, AD.fho_e_m, ...
            'n_steps', n_steps);
        s_next = stat_weight_mode(q + 1, mode);
        k_up = k_down_next * (s_next / s_q) * exp(-eps_m / (AD.k * T));
    end
    Z = Z + s_q * boltz;
    acc = acc + s_q * boltz * (k_up - k_down);
end
S = acc / Z;
end

function sw = stat_weight_mode(level, mode)
switch mode
    case 1
        sw = 1;
    case 2
        sw = level + 1;
    case {3, 4}
        sw = 0.5 * (level + 1) * (level + 2);
    otherwise
        error('main_test:bad_mode', 'mode must be 1, 2, 3, or 4.');
end
end

function r = safe_ratio(a, b)
if abs(b) <= realmin
    r = NaN;
else
    r = a / b;
end
end

function dy = rhs_simplified_2t(~, y, AD, processes, cfg)
processes = normalize_process_list(processes);
T = max(y(1) * AD.T0, 1e-6);
Tv = max(y(2) * AD.T0, 1e-6);

Rvibr = direct_vt_Rvibr_2t(T, Tv, AD, cfg, processes);

vv_processes = vv_only_processes(processes);
if ~isempty(vv_processes)
    terms = ch4_relax_terms_simplified(T, Tv, AD, ...
        'processes', vv_processes, 'n_steps_vv', cfg.n_steps_vv);
    Rvibr = Rvibr + terms.Rvibr;
end

A = capacity_matrix_2t(Tv, AD);
B = [0; Rvibr];
dy = A \ B;
end

function dy = rhs_simplified_3t(~, y, AD, processes, cfg)
processes = normalize_process_list(processes);
T = max(y(1) * AD.T0, 1e-6);
T13 = max(y(2) * AD.T0, 1e-6);
T24 = max(y(3) * AD.T0, 1e-6);

R13 = 0;
R24 = direct_vt_R24_3t(T, T24, AD, cfg, processes);

vv_processes = vv_only_processes(processes);
if ~isempty(vv_processes)
    terms = ch4_relax_terms_simplified(T, T13, T24, AD, ...
        'processes', vv_processes, 'n_steps_vv', cfg.n_steps_vv);
    R13 = R13 + terms.R13;
    R24 = R24 + terms.R24;
end

A = capacity_matrix_3t(T13, T24, AD);
B = [0; R13; R24];
dy = A \ B;
end

function Rvibr = direct_vt_Rvibr_2t(T, Tv, AD, cfg, processes)
eps = mode_energies(AD);
scale_nd = AD.n0 * AD.tau / (AD.k * AD.T0);
Rvibr = 0;
if has_named_process(processes, 'VT2')
    S2 = direct_vt_sum_mode(T, Tv, 2, eps.eps2, AD, cfg.n_steps_vt);
    Rvibr = Rvibr + scale_nd * eps.eps2 * S2;
end
if has_named_process(processes, 'VT4')
    S4 = direct_vt_sum_mode(T, Tv, 4, eps.eps4, AD, cfg.n_steps_vt);
    Rvibr = Rvibr + scale_nd * eps.eps4 * S4;
end
end

function R24 = direct_vt_R24_3t(T, T24, AD, cfg, processes)
eps = mode_energies(AD);
scale_nd = AD.n0 * AD.tau / (AD.k * AD.T0);
R24 = 0;
if has_named_process(processes, 'VT2')
    S2 = direct_vt_sum_mode(T, T24, 2, eps.eps2, AD, cfg.n_steps_vt);
    R24 = R24 + scale_nd * eps.eps2 * S2;
end
if has_named_process(processes, 'VT4')
    S4 = direct_vt_sum_mode(T, T24, 4, eps.eps4, AD, cfg.n_steps_vt);
    R24 = R24 + scale_nd * eps.eps4 * S4;
end
end

function processes_vv = vv_only_processes(processes)
processes = normalize_process_list(processes);
processes_vv = {};
for i = 1:numel(processes)
    name = char(string(processes{i}));
    if ~strcmpi(name, 'VT2') && ~strcmpi(name, 'VT4')
        processes_vv{end + 1} = name; %#ok<AGROW>
    end
end
end

function tf = has_named_process(processes, target)
processes = normalize_process_list(processes);
tf = false;
for i = 1:numel(processes)
    if strcmpi(char(string(processes{i})), target)
        tf = true;
        return
    end
end
end

function processes = normalize_process_list(processes)
if isempty(processes)
    processes = {};
    return
end
if ischar(processes) || isstring(processes)
    processes = cellstr(string(processes));
    return
end
while iscell(processes) && isscalar(processes) && iscell(processes{1})
    processes = processes{1};
end
if ~iscell(processes)
    processes = {processes};
end
end

function A = capacity_matrix_2t(Tv, AD)
w = AD.stw(:);
E = AD.e1234(:);
xi = -E ./ (AD.k * Tv);
x = normalized_boltzmann(w, xi);
q = E ./ (AD.k * Tv);
cv = sum(x .* q.^2) - sum(x .* q)^2;

A = eye(2);
A(1, 1) = 3;
A(1, 2) = cv;
A(2, 2) = cv;
end

function A = capacity_matrix_3t(T13, T24, AD)
eps = mode_energies(AD);
I0 = AD.inds(:, 1);
J0 = AD.inds(:, 2);
K0 = AD.inds(:, 3);
L0 = AD.inds(:, 4);

E13 = I0 * eps.eps1 + K0 * eps.eps3;
E24 = J0 * eps.eps2 + L0 * eps.eps4;

xi = -(E13 ./ (AD.k * T13) + E24 ./ (AD.k * T24));
x = normalized_boltzmann(AD.stw(:), xi);

E13m = sum(x .* E13);
E24m = sum(x .* E24);
varE13 = sum(x .* E13.^2) - E13m^2;
varE24 = sum(x .* E24.^2) - E24m^2;
covE13E24 = sum(x .* E13 .* E24) - E13m * E24m;

de13_dT13 = varE13 / (AD.k^2 * T13^2);
de13_dT24 = covE13E24 / (AD.k^2 * T24^2);
de24_dT13 = covE13E24 / (AD.k^2 * T13^2);
de24_dT24 = varE24 / (AD.k^2 * T24^2);

A = zeros(3, 3);
A(1, 1) = 3;
A(1, 2) = de13_dT13 + de24_dT13;
A(1, 3) = de13_dT24 + de24_dT24;
A(2, 2) = de13_dT13;
A(2, 3) = de13_dT24;
A(3, 2) = de24_dT13;
A(3, 3) = de24_dT24;
end

function x = normalized_boltzmann(w, xi)
shift = max(xi);
f = w .* exp(xi - shift);
x = f ./ sum(f);
end

function eps = mode_energies(AD)
eps.eps1 = AD.e1000 - AD.e0000;
eps.eps2 = AD.e0100 - AD.e0000;
eps.eps3 = AD.e0010 - AD.e0000;
eps.eps4 = AD.e0001 - AD.e0000;
end

function result = pack_result(run_case, X, Y, AD, elapsed)
result = empty_result();
result.name = run_case.name;
result.kind = run_case.kind;
result.dim = run_case.dim;
result.processes = run_case.processes;
result.time = X * AD.tau;
result.T = Y(:, 1) * AD.T0;
result.elapsed = elapsed;

if run_case.dim == 3
    result.T13 = Y(:, 2) * AD.T0;
    result.T24 = Y(:, 3) * AD.T0;
else
    result.Tv = Y(:, 2) * AD.T0;
end
end

function result = empty_result()
result = struct( ...
    'name', '', ...
    'kind', '', ...
    'dim', [], ...
    'processes', {{}}, ...
    'time', [], ...
    'T', [], ...
    'Tv', [], ...
    'T13', [], ...
    'T24', [], ...
    'elapsed', []);
end

function print_result_summary(result, T_eq)
fprintf('  elapsed: %.2f s\n', result.elapsed);
fprintf('  final T: %.6f K\n', result.T(end));
if result.dim == 3
    fprintf('  final T13: %.6f K, final T24: %.6f K\n', result.T13(end), result.T24(end));
else
    fprintf('  final Tv: %.6f K\n', result.Tv(end));
end
fprintf('  max final |T_i - T_eq|: %.3e K\n', equilibrium_error(result, T_eq));
end

function print_equilibrium_table(results, T_eq)
fprintf('\nFinal equilibrium check, T_eq = %.6f K\n', T_eq);
fprintf('%3s  %-38s  %10s  %10s  %10s  %12s\n', ...
    '#', 'case', 'T', 'Tv/T13', 'T24', 'max err K');
fprintf('%s\n', repmat('-', 1, 92));
for i = 1:numel(results)
    r = results(i);
    if r.dim == 3
        fprintf('%3d  %-38s  %10.4f  %10.4f  %10.4f  %12.3e\n', ...
            i, r.name, r.T(end), r.T13(end), r.T24(end), equilibrium_error(r, T_eq));
    else
        fprintf('%3d  %-38s  %10.4f  %10.4f  %10s  %12.3e\n', ...
            i, r.name, r.T(end), r.Tv(end), '-', equilibrium_error(r, T_eq));
    end
end
end

function err = equilibrium_error(result, T_eq)
temps = result.T(end);
if result.dim == 3
    temps = [temps, result.T13(end), result.T24(end)];
else
    temps = [temps, result.Tv(end)];
end
err = max(abs(temps - T_eq));
end

function plot_main_style(results, T_eq, cfg)
figure('Name', 'Main comparison');
hold on;

idx = 1:numel(results);
style = main_plot_styles();
legend_entries = {};

for i = 1:numel(idx)
    r = results(idx(i));
    st = style(i);
    [t, T_plot] = plotted_series(r.time, r.T);
    plot_temperature(t, T_plot, st.TLine, st.TColor, 'none');
    legend_entries{end + 1} = ['T - ' r.name]; %#ok<AGROW>

    if r.dim == 3
        [~, T13_plot] = plotted_series(r.time, r.T13);
        [~, T24_plot] = plotted_series(r.time, r.T24);
        plot_temperature(t, T13_plot, '--', st.T13Color, 'none');
        plot_temperature(t, T24_plot, '-', st.T24Color, 'none');
        legend_entries{end + 1} = ['T13 - ' r.name]; %#ok<AGROW>
        legend_entries{end + 1} = ['T24 - ' r.name]; %#ok<AGROW>
    else
        [~, Tv_plot] = plotted_series(r.time, r.Tv);
        plot_temperature(t, Tv_plot, '--', st.TvColor, 'none');
        legend_entries{end + 1} = ['Tv - ' r.name]; %#ok<AGROW>
    end
end

yline(T_eq, ':', sprintf('T_{eq}=%.2f K', T_eq), 'Color', [0.25, 0.25, 0.25], 'LineWidth', 1.2);
xline(cfg.plot_t_min, ':', 'start', 'Color', [0.35, 0.35, 0.35]);
xlim([cfg.plot_t_min, cfg.t_fin]);
ylim([min(250, cfg.Tv0 - 40), max(cfg.T0 + 40, T_eq + 40)]);
xlabel('t [s]');
ylabel('Temperature [K]');
title(sprintf('T_0 = %.0f K, T_{v0} = %.0f K', cfg.T0, cfg.Tv0));
legend(legend_entries, 'Location', 'northeast');
grid on;
set(gca, 'XScale', 'log');
set(gca, 'XMinorGrid', 'on');
end

function styles = main_plot_styles()
styles = struct( ...
    'TColor', { ...
        [0.00, 0.00, 0.00], ...
        [0.05, 0.55, 0.05], ...
        [0.55, 0.00, 0.15], ...
        [0.38, 0.10, 0.65], ...
        [0.00, 0.45, 0.55]}, ...
    'TvColor', { ...
        [0.00, 0.00, 0.00], ...
        [0.25, 0.65, 0.85], ...
        [0.25, 0.65, 0.85], ...
        [0.25, 0.65, 0.85], ...
        [0.25, 0.65, 0.85]}, ...
    'T13Color', { ...
        [0.00, 0.25, 0.85], ...
        [0.00, 0.25, 0.85], ...
        [0.00, 0.25, 0.85], ...
        [0.45, 0.10, 0.75], ...
        [0.10, 0.35, 0.75]}, ...
    'T24Color', { ...
        [0.85, 0.30, 0.05], ...
        [0.85, 0.30, 0.05], ...
        [0.85, 0.30, 0.05], ...
        [0.00, 0.55, 0.65], ...
        [0.90, 0.45, 0.00]}, ...
    'TLine', {'-', '-', '-', '-', '-'} ...
);
end

function plot_split_style(results, T_eq, cfg, title_text)
idx = 1:numel(results);
figure('Name', [title_text ' split']);
tiledlayout(2, 1, 'TileSpacing', 'compact');

nexttile;
hold on;
colors = lines(numel(idx));
legend_entries = {};
for i = 1:numel(idx)
    r = results(idx(i));
    [t, T_plot] = plotted_series(r.time, r.T);
    semilogx(t, T_plot, 'LineStyle', '-', 'Color', colors(i, :), 'LineWidth', 1.6);
    legend_entries{end + 1} = r.name; %#ok<AGROW>
end
yline(T_eq, 'k:', 'LineWidth', 1.1);
xlim([cfg.plot_t_min, cfg.t_fin]);
ylabel('T [K]');
title('Translational temperature');
legend(legend_entries, 'Location', 'best');
grid on;
set(gca, 'XMinorGrid', 'on');

nexttile;
hold on;
legend_entries = {};
for i = 1:numel(idx)
    r = results(idx(i));
    [t, ~] = plotted_series(r.time, r.T);
    if r.dim == 3
        [~, T13_plot] = plotted_series(r.time, r.T13);
        [~, T24_plot] = plotted_series(r.time, r.T24);
        semilogx(t, T13_plot, 'LineStyle', '--', 'Color', colors(i, :), 'LineWidth', 1.4);
        semilogx(t, T24_plot, 'LineStyle', '-.', 'Color', colors(i, :), 'LineWidth', 1.4);
        legend_entries{end + 1} = ['T13 - ' r.name]; %#ok<AGROW>
        legend_entries{end + 1} = ['T24 - ' r.name]; %#ok<AGROW>
    else
        [~, Tv_plot] = plotted_series(r.time, r.Tv);
        semilogx(t, Tv_plot, 'LineStyle', '--', 'Color', colors(i, :), 'LineWidth', 1.4);
        legend_entries{end + 1} = ['Tv - ' r.name]; %#ok<AGROW>
    end
end
yline(T_eq, 'k:', 'LineWidth', 1.1);
xlim([cfg.plot_t_min, cfg.t_fin]);
xlabel('t [s]');
ylabel('Vibrational temperatures [K]');
title(title_text);
legend(legend_entries, 'Location', 'best');
grid on;
set(gca, 'XMinorGrid', 'on');
end

function plot_temperature(t, y, line_style, color, marker)
semilogx(t, y, 'LineStyle', line_style, 'Color', color, 'LineWidth', 1.7);
if ~strcmp(marker, 'none')
    semilogx(t(1), y(1), marker, ...
        'Color', color, 'MarkerFaceColor', 'w', 'MarkerSize', 6, 'LineWidth', 1.2, ...
        'HandleVisibility', 'off');
    semilogx(t(end), y(end), marker, ...
        'Color', color, 'MarkerFaceColor', color, 'MarkerSize', 6, 'LineWidth', 1.2, ...
        'HandleVisibility', 'off');
end
end

function [t_plot, y_plot] = plotted_series(t, y)
mask = t(:) > 0;
t_plot = t(mask);
y_plot = y(mask);
end

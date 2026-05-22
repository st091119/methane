clear

% Compare 3T STS relaxation: all processes vs without VV32d and VV32-2d.

addpath(fullfile(fileparts(mfilename('fullpath')), 'fho_model'));

%% Settings
p0 = 101325;
T0 = 1000;
T13_0 = 700;
T24_0 = 300;
t_fin = 1;

options = odeset('RelTol', 1e-8, 'AbsTol', 1e-8);

%% Model data
AD = input_data;
fprintf('Spectrum: %s, lch4 = [%s], klop = %d\n', AD.spectrum_mode, ...
    num2str(AD.lch4), AD.klop);

n0 = p0 / (AD.k * T0);
sigma0 = pi * AD.r0^2;
tau = (4 * n0 * sigma0 * sqrt(AD.k * T0 / (pi * AD.m)))^(-1);
fprintf('Mean time between collisions tau = %.6e sec\n', tau);

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
AD.use_vv34_4d = true;
AD.use_vv34_2d = true;
AD.use_vv32_4d = true;

AD = ch4_relax_topology(AD);

tspan = [0, t_fin] ./ tau;
Y0_3t = [1; T13_0 / T0; T24_0 / T0];

cases = struct( ...
    'name', {'full', 'no32'}, ...
    'label', {'all processes', 'no VV32d, VV32-2d'}, ...
    'use_vv32d', {true, false}, ...
    'use_vv32_2d', {true, false} ...
);

results = struct();

for ic = 1:numel(cases)
    AD_run = AD;
    AD_run.use_vv32d = cases(ic).use_vv32d;
    AD_run.use_vv32_2d = cases(ic).use_vv32_2d;

    fprintf('\n--- Solving 3T STS: %s (%s) ---\n', cases(ic).name, cases(ic).label);
    tic;
    [X, Y] = ode15s(@(t, y) rpart_mt_3t_sts(t, y, AD_run), tspan, Y0_3t, options);
    elapsed = toc;
    fprintf('  Done in %.2f sec, %d time points\n', elapsed, numel(X));

    results(ic).elapsed = elapsed;
    results(ic).time = X * tau;
    results(ic).T = Y(:, 1) * T0;
    results(ic).T13 = Y(:, 2) * T0;
    results(ic).T24 = Y(:, 3) * T0;
    results(ic).label = cases(ic).label;

    fprintf('  Final T  = %.4f K\n', results(ic).T(end));
    fprintf('  Final T13 = %.4f K\n', results(ic).T13(end));
    fprintf('  Final T24 = %.4f K\n', results(ic).T24(end));
end

t_full = results(1).elapsed;
t_no32 = results(2).elapsed;
dt_solve = t_full - t_no32;

fprintf('\n=== Время расчёта (ode15s) ===\n');
fprintf('  Полная схема (все процессы):     %.2f sec\n', t_full);
fprintf('  Упрощённая (без VV32d, VV32-2d): %.2f sec\n', t_no32);
fprintf('  Разница (полная - упрощённая):   %.2f sec', dt_solve);
if t_no32 > 0
    fprintf(' (ускорение упрощённой x%.2f)\n', t_full / t_no32);
else
    fprintf('\n');
end

%% Common time grid and differences
t_common = union(results(1).time, results(2).time);

T_full = interp1(results(1).time, results(1).T, t_common);
T_no32 = interp1(results(2).time, results(2).T, t_common);
T13_full = interp1(results(1).time, results(1).T13, t_common);
T13_no32 = interp1(results(2).time, results(2).T13, t_common);
T24_full = interp1(results(1).time, results(1).T24, t_common);
T24_no32 = interp1(results(2).time, results(2).T24, t_common);

dT = T_full - T_no32;
dT13 = T13_full - T13_no32;
dT24 = T24_full - T24_no32;

%% Metrics
names = {'T', 'T13', 'T24'};
d_all = {dT, dT13, dT24};
ref_end = [T_full(end), T13_full(end), T24_full(end)];

fprintf('\n=== Difference metrics (full - no VV32) ===\n');
fprintf('%-8s %12s %12s %12s %12s\n', 'var', 'max_abs[K]', 'rms[K]', 'final[K]', 'final_rel[%]');
max_over_all = 0;
for iv = 1:3
    dv = d_all{iv};
    max_abs = max(abs(dv));
    rms_val = sqrt(mean(dv.^2));
    final_abs = abs(dv(end));
    final_rel = 100 * final_abs / max(ref_end(iv), 1);
    max_over_all = max(max_over_all, max_abs);
    fprintf('%-8s %12.6f %12.6f %12.6f %12.6f\n', ...
        names{iv}, max_abs, rms_val, final_abs, final_rel);
end
fprintf('max_over_all (worst |Delta| over T, T13, T24) = %.6f K\n', max_over_all);

%% Figure 1: temperature trajectories
figure('Name', '3T compare: temperatures');
hold on;
plot(results(1).time, results(1).T, 'k-', 'LineWidth', 2);
plot(results(2).time, results(2).T, 'k--', 'LineWidth', 2);
plot(results(1).time, results(1).T13, 'r-', 'LineWidth', 2);
plot(results(2).time, results(2).T13, 'r--', 'LineWidth', 2);
plot(results(1).time, results(1).T24, 'b-', 'LineWidth', 2);
plot(results(2).time, results(2).T24, 'b--', 'LineWidth', 2);
set(gca, 'XScale', 'log');
xlabel('t [sec]');
ylabel('Temperature [K]');
legend({'T — full', 'T — no VV32', 'T_{13} — full', 'T_{13} — no VV32', ...
    'T_{24} — full', 'T_{24} — no VV32'}, 'Location', 'best');
title(sprintf('3T STS: T_0=%g K, T_{13,0}=%g K, T_{24,0}=%g K', T0, T13_0, T24_0));
grid on;
box on;

%% Figure 2: absolute and relative differences
figure('Name', '3T compare: differences');
subplot(3, 1, 1);
semilogx(t_common, abs(dT), 'k-', 'LineWidth', 2);
ylabel('|\\DeltaT| [K]');
title('|\\Delta| = full - no VV32');
grid on;

subplot(3, 1, 2);
semilogx(t_common, abs(dT13), 'r-', 'LineWidth', 2);
ylabel('|\\DeltaT_{13}| [K]');
grid on;

subplot(3, 1, 3);
semilogx(t_common, abs(dT24), 'b-', 'LineWidth', 2);
ylabel('|\\DeltaT_{24}| [K]');
xlabel('t [sec]');
grid on;

figure('Name', '3T compare: relative differences');
subplot(3, 1, 1);
semilogx(t_common, 100 * abs(dT) ./ max(T_full, 1), 'k-', 'LineWidth', 2);
ylabel('rel. |\\DeltaT| [%]');
title('Relative difference vs full case');
grid on;

subplot(3, 1, 2);
semilogx(t_common, 100 * abs(dT13) ./ max(T13_full, 1), 'r-', 'LineWidth', 2);
ylabel('rel. |\\DeltaT_{13}| [%]');
grid on;

subplot(3, 1, 3);
semilogx(t_common, 100 * abs(dT24) ./ max(T24_full, 1), 'b-', 'LineWidth', 2);
ylabel('rel. |\\DeltaT_{24}| [%]');
xlabel('t [sec]');
grid on;

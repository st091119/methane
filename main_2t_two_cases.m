clear

addpath(fullfile(fileparts(mfilename('fullpath')), 'fho_model'));

%% Two-temperature test cases
p0 = 101325;
t_fin = 1;
options = odeset('RelTol', 1e-8, 'AbsTol', 1e-8);

test_cases = struct( ...
    'name', {'TC1', 'TC2'}, ...
    'T0',   {1000, 300}, ...
    'Tv0',  {300, 900} ...
);

run_models = struct( ...
    'name',  {'Wang-Springer', 'Landau-Teller (FHO)', 'Hybrid 2T'}, ...
    'rp',    {'rpart_mt_lt',   'rpart_mt_lt',         'rpart_mt_sts'}, ...
    'sw_rt', {'vt_rel_time_wang', 'fho',              'fho'} ...
);

for ic = 1:numel(test_cases)
    T0 = test_cases(ic).T0;
    Tv0 = test_cases(ic).Tv0;

    AD = setup_case_ad(T0, p0);
    tspan = [0, t_fin] ./ AD.tau;
    Y0 = [1; Tv0 / T0];

    fprintf('\n=== %s: T0 = %.0f K, Tv0 = %.0f K ===\n', ...
        test_cases(ic).name, T0, Tv0);

    results = struct();
    for im = 1:numel(run_models)
        AD_run = AD;
        AD_run.sw_rt = run_models(im).sw_rt;
        RP = str2func(run_models(im).rp);

        fprintf('--- Solving %s ---\n', run_models(im).name);
        [X, Y] = ode15s(@(t, y) RP(t, y, AD_run), tspan, Y0, options);

        results(im).name = run_models(im).name;
        results(im).time = X * AD.tau;
        results(im).T = Y(:, 1) * T0;
        results(im).Tv = Y(:, 2) * T0;

        fprintf('  Final T  = %.6g K\n', results(im).T(end));
        fprintf('  Final Tv = %.6g K\n', results(im).Tv(end));
    end

    plot_case_results(results, T0, Tv0, test_cases(ic).name);
end

function AD = setup_case_ad(T0, p0)
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

AD.fho_alpha = 5.174e10;
AD.fho_e_m = 776.4;
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
end

function plot_case_results(results, T0, Tv0, case_name)
colors = lines(numel(results));

figure;
hold on;
legend_entries = {};
for im = 1:numel(results)
    semilogx(results(im).time, results(im).T, ...
        'Color', colors(im, :), 'LineStyle', '-', 'LineWidth', 2);
    semilogx(results(im).time, results(im).Tv, ...
        'Color', colors(im, :), 'LineStyle', '--', 'LineWidth', 2);

    legend_entries{end + 1} = ['T - ' results(im).name]; %#ok<AGROW>
    legend_entries{end + 1} = ['T_v - ' results(im).name]; %#ok<AGROW>
end

set(gca, 'XScale', 'log');
xlabel('t [sec]');
ylabel('Temperature [K]');
title(sprintf('%s: T_0 = %.0f K, T_{v0} = %.0f K', case_name, T0, Tv0));
legend(legend_entries, 'Location', 'best');
grid on;
end

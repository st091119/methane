% test_fho.m — MATLAB test for FHO model: VT rates and relaxation times
% Equivalent to the Python reference test using Fit C kinetic parameters
% Run this script from the fho_model/ folder or from the parent folder.

clear; clc;

% Add paths
thisDir   = fileparts(mfilename('fullpath'));
parentDir = fullfile(thisDir, '..');
addpath(thisDir);
addpath(parentDir);

%% ══════════════════════════════════════════════════════════════
%  FHO model parameters, obtained by fitting experimental data for CH4-CH4 collisions
%  ══════════════════════════════════════════════════════════════
alpha_test   = 4.6947967158e+10;   % [m^-1] - Morse range parameter
e_m_test     = 9.5648964013e+02;   % [K]    - Morse well depth / k_B
steric2_test = 6.4499971926e-02;   % steric factor for mode 2
steric4_test = 6.7853797863e-03;   % steric factor for mode 4

fprintf('Parameters for MATLAB test:\n');
fprintf('  alpha   = %.10e m^-1\n', alpha_test);
fprintf('  e_m     = %.10e K\n',    e_m_test);
fprintf('  steric2 = %.10e\n',      steric2_test);
fprintf('  steric4 = %.10e\n',      steric4_test);
fprintf('\n');

%% ── 1. VT rate constants: 4 test cases ──────────────────────
% Each row: {T, qi, qf, mode (1-based), steric, label}
test_cases = {
    300,  1, 0, 2, steric2_test, 'VT2: (0,1,0,0)->(0,0,0,0) T=300K';
    1000, 2, 1, 2, steric2_test, 'VT2: (0,2,0,0)->(0,1,0,0) T=1000K';
    300,  1, 0, 4, steric4_test, 'VT4: (0,0,0,1)->(0,0,0,0) T=300K';
    1000, 3, 2, 4, steric4_test, 'VT4: (0,0,0,3)->(0,0,0,2) T=1000K';
};

n_cases = size(test_cases, 1);

fprintf('%s\n', repmat('=', 1, 75));
fprintf('%-45s %14s  %14s\n', 'Case', 'k_VT [m3/s]', 'k_VT [cm3/s]');
fprintf('%s\n', repmat('-', 1, 75));

for ic = 1:n_cases
    T_c    = test_cases{ic, 1};
    qi_c   = test_cases{ic, 2};
    qf_c   = test_cases{ic, 3};
    mode_c = test_cases{ic, 4};
    st_c   = test_cases{ic, 5};
    label  = test_cases{ic, 6};

    k_val = rates_fho_vt(T_c, qi_c, qf_c, mode_c, st_c, ...
                         alpha_test, e_m_test, ...
                         'g0_max', 80.0, 'n_steps', 10000);

    fprintf('  %-43s %14.6e  %14.6e\n', label, k_val, k_val * 1e6);
end

%% ── 2. Relaxation times: p*tau at selected temperatures ─────
T_list = [300, 500, 1000, 2000];

fprintf('\n');
fprintf('%s\n', repmat('=', 1, 75));
fprintf('%8s  %14s  %14s  %16s\n', 'T [K]', 'ptau2 [atm*s]', 'ptau4 [atm*s]', 'ptau_tot [atm*s]');
fprintf('%s\n', repmat('-', 1, 75));

for iT = 1:length(T_list)
    T_r = T_list(iT);

    res = relaxation_time_kinetic(T_r, alpha_test, e_m_test, ...
                                  steric2_test, steric4_test, ...
                                  'n_steps', 10000);

    fprintf('  %6d  %14.6e  %14.6e  %16.6e\n', ...
            T_r, res.ptau2, res.ptau4, res.ptau_total);
end

%% Print relative populations at T range from T_list
AD = input_data();
fprintf('\nRelative populations of mode 2 and mode 4 at selected temperatures:\n');
fprintf('%8s  %14s\n', 'T [K]', 'alpha');
for iT = 1:length(T_list)
    T_r = T_list(iT);
    pop = relative_boltzmann_population_full(T_r, AD);
    fprintf('  %6d  %14.6e\n', T_r, pop);
end

fprintf('\nDone.\n');

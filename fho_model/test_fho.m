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
alpha2_test   = 4.3479e10;   % [m^-1] - Morse range parameter
e_m2_test     = 1434.2;      % [K]    - Morse well depth / k_B
alpha4_test   = 6.6145e10;   % [m^-1] - Morse range parameter
e_m4_test     = 280.0;       % [K]    - Morse well depth / k_B
steric2_test  = 0.0397;      % steric factor for mode 2
steric4_test  = 0.0050;      % steric factor for mode 4

fprintf('Parameters for MATLAB test:\n');
fprintf('  alpha 2nd mode   = %.10e m^-1\n', alpha2_test);
fprintf('  e_m 2nd mode     = %.10e K\n',    e_m2_test);
fprintf('  alpha 4th mode   = %.10e m^-1\n', alpha4_test);
fprintf('  e_m 4th mode     = %.10e K\n',    e_m4_test);
fprintf('  steric 2nd mode  = %.10e\n',      steric2_test);
fprintf('  steric 4th mode  = %.10e\n',      steric4_test);
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
    
    if mode_c == 2
        alpha_c = alpha2_test;
        e_m_c   = e_m2_test;
    elseif mode_c == 4
        alpha_c = alpha4_test;
        e_m_c   = e_m4_test;
    else
        error('Invalid mode in test case: %d', mode_c);
    end
    k_val = rates_fho_vt(T_c, qi_c, qf_c, mode_c, st_c, ...
                         alpha_c, e_m_c, ...
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

    res = relaxation_time_kinetic(T_r, alpha2_test, e_m2_test, ...
                                  alpha4_test, e_m4_test, ...
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

%% ══════════════════════════════════════════════════════════════
%  VV rate coefficients — 3 processes (ν₃→ν₄ deactivation)
%  ══════════════════════════════════════════════════════════════
%  Morse parameters: ν₄ mode from Fit C
%  Steric factors (VT and VV) obtained from fitting

% ── Steric factors ──
steric_vt_v3   = sqrt(steric2_test * steric4_test);  % ≈ 0.0209 (geometric mean)
steric_vt_v3_4 = 0.99;                               % near-resonant ν₃→2ν₄
steric_vv_34s   = 0.068;    % VV_{3-4}^s   (intramolecular ν₃→ν₄)
steric_vv_34d   = 1.0;      % VV_{3-4}^d   (intramolecular ν₃→2ν₄)
steric_vv_34_4d = 0.2;      % VV_{3-4,4}^d (intermolecular ν₃→ν₄+ν₄')

fprintf('\n');
fprintf('VV steric factors:\n');
fprintf('  steric_vt_v3   = %.6e  (mean of steric2 and steric4 for convenience)\n', steric_vt_v3);
fprintf('  steric_vt_v3_4 = %.2f\n', steric_vt_v3_4);
fprintf('  steric_vv_34s  = %.3f\n', steric_vv_34s);
fprintf('  steric_vv_34d  = %.3f\n', steric_vv_34d);
fprintf('  steric_vv_34_4d= %.3f\n', steric_vv_34_4d);

%  Three VV processes:
%   VV_{3-4}^s, example   (0,0,1,0)->(0,0,0,1), partner (0,2,0,1)->(0,2,0,1)  [intramolecular]
%   VV_{3-4}^d, example   (0,0,1,0)->(0,0,0,2), partner (0,0,0,0)->(0,0,0,0)  [intramolecular]
%   VV_{3-4,4}^d, example (0,0,1,0)->(0,0,0,1), partner (0,0,0,0)->(0,0,0,1)  [intermolecular]

vv_cases = {
    % state_i        state_f        state_k        state_kf       svt             svv             label
    [0 0 1 0],  [0 0 0 1],  [0 2 0 1],  [0 2 0 1],  steric_vt_v3,   steric_vv_34s,   'VV34s';
    [0 0 1 0],  [0 0 0 2],  [0 0 0 0],  [0 0 0 0],  steric_vt_v3_4, steric_vv_34d,   'VV34d';
    [0 0 1 0],  [0 0 0 1],  [0 0 0 0],  [0 0 0 1],  steric_vt_v3,   steric_vv_34_4d, 'VV34_4d';
};

%% ── 3. VV rate coefficients at selected temperatures ────────
T_vv_list = [300, 500, 1000, 1400];

fprintf('\n');
fprintf('%s\n', repmat('=', 1, 85));
fprintf('%-12s %6s  %14s  %14s\n', 'Process', 'T [K]', 'k_VV [m3/s]', 'k_VV [cm3/s]');
fprintf('%s\n', repmat('-', 1, 85));

for ic = 1:size(vv_cases, 1)
    si  = vv_cases{ic, 1};
    sf  = vv_cases{ic, 2};
    sk  = vv_cases{ic, 3};
    skf = vv_cases{ic, 4};
    svt = vv_cases{ic, 5};
    svv = vv_cases{ic, 6};
    lbl = vv_cases{ic, 7};

    for iT = 1:length(T_vv_list)
        T_v = T_vv_list(iT);
        k_val = rates_fho_vv(T_v, si, sf, sk, skf, svt, svv, ...
                             alpha4_test, e_m4_test, ...
                             'g0_max', 80.0, 'n_steps', 10000);
        fprintf('  %-10s %6d  %14.6e  %14.6e\n', lbl, T_v, k_val, k_val * 1e6);
    end
    fprintf('\n');
end

%% ── 4. VV relaxation times p*tau (simplified kT/k10) ───────
fprintf('%s\n', repmat('=', 1, 85));
fprintf('%-12s %6s  %14s  %14s\n', 'Process', 'T [K]', 'k10 [m3/s]', 'ptau [atm*s]');
fprintf('%s\n', repmat('-', 1, 85));

k_B_val = AD.k;

for ic = 1:size(vv_cases, 1)
    si  = vv_cases{ic, 1};
    sf  = vv_cases{ic, 2};
    sk  = vv_cases{ic, 3};
    skf = vv_cases{ic, 4};
    svt = vv_cases{ic, 5};
    svv = vv_cases{ic, 6};
    lbl = vv_cases{ic, 7};

    for iT = 1:length(T_vv_list)
        T_v = T_vv_list(iT);

        % k_{1->0} rate (donor mode changes 1->0 in mode 3)
        k10 = rates_fho_vv(T_v, si, sf, sk, skf, svt, svv, ...
                           alpha4_test, e_m4_test, ...
                           'g0_max', 80.0, 'n_steps', 10000);

        if k10 > 0
            ptau_val = k_B_val * T_v / k10 / 101325;
        else
            ptau_val = Inf;
        end
        fprintf('  %-10s %6d  %14.6e  %14.6e\n', lbl, T_v, k10, ptau_val);
    end
    fprintf('\n');
end

%% ── 5. VV relaxation times — kinetic theory (integrated) ────
%  Uses relaxation_time_kinetic with 'vv_cases' parameter.
%  vv_cases format: {si0, sf0, sk, skf, svt, svv, alpha_vv, e_m_vv, label}

vv_cases_kin = {
    [0 0 1 0], [0 0 0 1], [0 2 0 1], [0 2 0 1], steric_vt_v3,   steric_vv_34s,   alpha4_test, e_m4_test, 'VV34s';
    [0 0 1 0], [0 0 0 2], [0 0 0 0], [0 0 0 0], steric_vt_v3_4, steric_vv_34d,   alpha4_test, e_m4_test, 'VV34d';
    [0 0 1 0], [0 0 0 1], [0 0 0 0], [0 0 0 1], steric_vt_v3,   steric_vv_34_4d, alpha4_test, e_m4_test, 'VV34_4d';
};

fprintf('\n');
fprintf('%s\n', repmat('=', 1, 100));
fprintf('%-12s %6s  %18s  %18s\n', 'Process', 'T [K]', 'ptau_kin [atm*s]', 'ptau_simpl [atm*s]');
fprintf('%s\n', repmat('-', 1, 100));

for iT = 1:length(T_vv_list)
    T_v = T_vv_list(iT);

    res_vv = relaxation_time_kinetic(T_v, alpha2_test, e_m2_test, ...
                                     alpha4_test, e_m4_test, ...
                                     steric2_test, steric4_test, ...
                                     'n_steps', 10000, ...
                                     'vv_cases', vv_cases_kin);

    for ic = 1:size(vv_cases_kin, 1)
        lbl = vv_cases_kin{ic, 9};
        fprintf('  %-10s %6d  %18.6e  %18.6e\n', lbl, T_v, ...
                res_vv.ptau_vv{ic}, res_vv.ptau_vv_simpl{ic});
    end
    fprintf('\n');
end

fprintf('\nDone.\n');

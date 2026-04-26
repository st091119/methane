% test_fho.m - MATLAB reference output for the CH4 FHO model.
% Run from fho_model/ or from the project root.

clear; clc;

thisDir = fileparts(mfilename('fullpath'));
parentDir = fullfile(thisDir, '..');
addpath(thisDir);
addpath(parentDir);

%% Fitted parameter set
alpha = 5.174e10;      % [m^-1]
e_m = 776.4;           % [K]
steric2_vt = 0.002612;
steric4_vt = 0.028546;

fprintf('Parameters for MATLAB FHO test:\n');
fprintf('  alpha      = %.10e m^-1\n', alpha);
fprintf('  e_m        = %.10e K\n', e_m);
fprintf('  steric VT2 = %.10e\n', steric2_vt);
fprintf('  steric VT4 = %.10e\n', steric4_vt);
fprintf('\n');

%% VT rate constants
vt_cases = {
    300,  1, 0, 2, steric2_vt, 'VT2: (0,1,0,0)->(0,0,0,0)';
    1000, 2, 1, 2, steric2_vt, 'VT2: (0,2,0,0)->(0,1,0,0)';
    300,  1, 0, 4, steric4_vt, 'VT4: (0,0,0,1)->(0,0,0,0)';
    1000, 3, 2, 4, steric4_vt, 'VT4: (0,0,0,3)->(0,0,0,2)';
};

fprintf('%s\n', repmat('=', 1, 90));
fprintf('VT rate constants\n');
fprintf('%s\n', repmat('-', 1, 90));
fprintf('%-44s %6s  %14s  %14s\n', 'Case', 'T [K]', 'k_VT [m3/s]', 'k_VT [cm3/s]');
fprintf('%s\n', repmat('-', 1, 90));

for ic = 1:size(vt_cases, 1)
    T_val = vt_cases{ic, 1};
    qi = vt_cases{ic, 2};
    qf = vt_cases{ic, 3};
    mode = vt_cases{ic, 4};
    steric = vt_cases{ic, 5};
    label = vt_cases{ic, 6};

    k_val = rates_fho_vt(T_val, qi, qf, mode, steric, alpha, e_m, ...
        'g0_max', 80.0, 'n_steps', 10000);

    fprintf('%-44s %6d  %14.6e  %14.6e\n', label, T_val, k_val, k_val * 1e6);
end

%% VT relaxation times
T_list = [300, 500, 1000, 2000];

fprintf('\n%s\n', repmat('=', 1, 90));
fprintf('VT relaxation references\n');
fprintf('%s\n', repmat('-', 1, 90));
fprintf('%8s  %14s  %14s  %16s\n', 'T [K]', 'ptau2 [atm s]', 'ptau4 [atm s]', 'ptau_tot [atm s]');
fprintf('%s\n', repmat('-', 1, 90));

for iT = 1:numel(T_list)
    T_val = T_list(iT);
    res = relaxation_time_kinetic(T_val, alpha, e_m, steric2_vt, steric4_vt, ...
                                  'n_steps', 10000);

    fprintf('%8d  %14.6e  %14.6e  %16.6e\n', ...
        T_val, res.ptau2, res.ptau4, res.ptau_total);
end

AD = input_data();
fprintf('\nRelative populations of mode 2 and mode 4:\n');
fprintf('%8s  %14s\n', 'T [K]', 'alpha_pop');
for iT = 1:numel(T_list)
    T_val = T_list(iT);
    fprintf('%8d  %14.6e\n', T_val, relative_boltzmann_population_full(T_val, AD));
end

%% Direct VV reactions and finalized steric factors
steric_vv_34s = 1.0;
steric_vv_34d = 0.06;
steric_vv_32d = 0.000625;
steric_vv_34_4d = 1.0;
steric_vv_44 = 0.419;
steric_vv_44_deact = 0.437;
steric_vv_24b = 0.99;
steric_vv_33 = 0.42;
steric_vv_13 = 0.99;

vv_reactions = {
    struct('key', '34s_vv',   'state_i', [0 0 1 0], 'state_f', [0 0 0 1], 'state_k', [0 0 0 0], 'state_kf', [0 0 0 0], 'steric_vv', steric_vv_34s,      'gamma', 0.5, 'params_name', 'params',       'formula', 'nu3 -> nu4',                    'note', 'cell 107 direct rate_vv');
    struct('key', '34d_vv',   'state_i', [0 0 1 0], 'state_f', [0 0 0 2], 'state_k', [0 0 0 0], 'state_kf', [0 0 0 0], 'steric_vv', steric_vv_34d,      'gamma', 1.0, 'params_name', 'params_vv34_d','formula', 'nu3 -> 2nu4',                   'note', 'cell 107 direct rate_vv');
    struct('key', '32d_vv',   'state_i', [0 0 1 0], 'state_f', [0 2 0 0], 'state_k', [0 0 0 0], 'state_kf', [0 0 0 0], 'steric_vv', steric_vv_32d,      'gamma', 0.5, 'params_name', 'params_vv32_d','formula', 'nu3 -> 2nu2',                   'note', 'cell 119 fitted');
    struct('key', '34_4d',    'state_i', [0 0 1 0], 'state_f', [0 0 0 1], 'state_k', [0 0 0 0], 'state_kf', [0 0 0 1], 'steric_vv', steric_vv_34_4d,    'gamma', 0.5, 'params_name', 'params',       'formula', 'nu3 -> nu4 + nu4_partner',      'note', 'cell 107 sharing channel');
    struct('key', '44_swap',  'state_i', [0 0 0 1], 'state_f', [0 0 0 0], 'state_k', [0 0 0 0], 'state_kf', [0 0 0 1], 'steric_vv', steric_vv_44,       'gamma', 0.5, 'params_name', 'params_vv44',  'formula', 'nu4 swap',                       'note', 'cell 111/113');
    struct('key', '44_deact', 'state_i', [0 0 0 2], 'state_f', [0 0 0 1], 'state_k', [0 0 0 0], 'state_kf', [0 0 0 1], 'steric_vv', steric_vv_44_deact, 'gamma', 0.5, 'params_name', 'params_vv44',  'formula', '2nu4 -> nu4 + nu4_partner',     'note', 'cell 111/113');
    struct('key', '24b',      'state_i', [0 1 0 0], 'state_f', [0 0 0 0], 'state_k', [0 0 0 0], 'state_kf', [0 0 0 1], 'steric_vv', steric_vv_24b,      'gamma', 0.5, 'params_name', 'params',       'formula', 'nu2 -> nu4_partner',            'note', 'cell 117');
    struct('key', '33_swap',  'state_i', [0 0 1 0], 'state_f', [0 0 0 0], 'state_k', [0 0 0 0], 'state_kf', [0 0 1 0], 'steric_vv', steric_vv_33,       'gamma', 0.5, 'params_name', 'params',       'formula', 'nu3 swap',                       'note', 'cell 119');
    struct('key', '13_swap',  'state_i', [1 0 0 0], 'state_f', [0 0 0 0], 'state_k', [0 0 0 0], 'state_kf', [0 0 1 0], 'steric_vv', steric_vv_13,       'gamma', 0.5, 'params_name', 'params',       'formula', 'nu1 -> nu3_partner',            'note', 'cell 119');
};

vv_lookup = containers.Map('KeyType', 'char', 'ValueType', 'any');
for ic = 1:numel(vv_reactions)
    vv_lookup(vv_reactions{ic}.key) = vv_reactions{ic};
end

fprintf('\n%s\n', repmat('=', 1, 150));
fprintf('SUMMARY OF DIRECT VV REACTIONS AND FINAL VV STERIC FACTORS\n');
fprintf('%s\n', repmat('=', 1, 150));
for ic = 1:numel(vv_reactions)
    rxn = vv_reactions{ic};
    fprintf('%-10s  CH4%s + CH4%s -> CH4%s + CH4%s    steric_vv = %.6f\n', ...
        rxn.key, state_str(rxn.state_i), state_str(rxn.state_k), ...
        state_str(rxn.state_f), state_str(rxn.state_kf), rxn.steric_vv);
end

%% Direct VV rates for Python/MATLAB comparison
T_vv_test = [300, 500, 1000, 1400];
matlab_rate_keys = {'34s_vv', '34d_vv', '32d_vv', '34_4d', ...
    '44_swap', '44_deact', '24b', '33_swap', '13_swap'};

fprintf('\n%s\n', repmat('=', 1, 120));
fprintf('Direct VV rate references for MATLAB (simplified)\n');
fprintf('%s\n', repmat('=', 1, 120));
fprintf('%-16s %6s  %14s  %-14s\n', 'Reaction', 'T [K]', 'k [m3/s]', 'params');
fprintf('%s\n', repmat('-', 1, 120));

for ik = 1:numel(matlab_rate_keys)
    rxn = vv_lookup(matlab_rate_keys{ik});
    for iT = 1:numel(T_vv_test)
        T_val = T_vv_test(iT);
        k_val = rates_fho_vv(T_val, rxn.state_i, rxn.state_f, rxn.state_k, rxn.state_kf, ...
            rxn.steric_vv, alpha, e_m, 'g0_max', 80.0, 'n_steps', 10000, ...
            'gamma', rxn.gamma);
        fprintf('%-16s %6d  %14.6e  %-14s\n', rxn.key, T_val, k_val, rxn.params_name);
    end
    fprintf('\n');
end

%% Mode-3 relaxation references
mode3_relaxation_keys = {'34s_vv', '34d_vv', '32d_vv', '34_4d'};

fprintf('%s\n', repmat('=', 1, 120));
fprintf('Mode-3 relaxation references used in p*tau checks\n');
fprintf('%s\n', repmat('=', 1, 120));
fprintf('%-12s %6s  %14s  %20s  %18s\n', ...
    'Reaction', 'T [K]', 'k_10 [m3/s]', 'ptau_simpl [atm s]', 'ptau_kin [atm s]');
fprintf('%s\n', repmat('-', 1, 120));

for ik = 1:numel(mode3_relaxation_keys)
    rxn = vv_lookup(mode3_relaxation_keys{ik});
    for iT = 1:numel(T_vv_test)
        T_val = T_vv_test(iT);
        [k10_val, ptau_simpl, ptau_kin] = compute_mode3_relaxation(rxn, T_val, alpha, e_m, AD);
        fprintf('%-12s %6d  %14.6e  %20.6e  %18.6e\n', ...
            rxn.key, T_val, k10_val, ptau_simpl, ptau_kin);
    end
    fprintf('\n');
end

fprintf('%s\n', repmat('=', 1, 120));
fprintf('Finalized VV steric factors used in the fitted set\n');
fprintf('%s\n', repmat('=', 1, 120));
fprintf('steric_vv_34s      = %.6f\n', steric_vv_34s);
fprintf('steric_vv_34d      = %.6f\n', steric_vv_34d);
fprintf('steric_vv_32d      = %.6f\n', steric_vv_32d);
fprintf('steric_vv_34_4d    = %.6f\n', steric_vv_34_4d);
fprintf('steric_vv_44       = %.6f\n', steric_vv_44);
fprintf('steric_vv_44_deact = %.6f\n', steric_vv_44_deact);
fprintf('steric_vv_24b      = %.6f\n', steric_vv_24b);
fprintf('steric_vv_33       = %.6f\n', steric_vv_33);
fprintf('steric_vv_13       = %.6f\n', steric_vv_13);
fprintf('Morse params       = alpha %.10e, e_m %.10e\n', alpha, e_m);

fprintf('\nDone.\n');

function text = state_str(state)
text = sprintf('(%d,%d,%d,%d)', state(1), state(2), state(3), state(4));
end

function [k10_val, ptau_simpl, ptau_kin] = compute_mode3_relaxation(rxn, T_val, alpha, e_m, AD)
ATM = 101325.0;
eps_3 = AD.h * AD.c * AD.omega(3);
max_level_3 = AD.lch4(3);

Z3 = 0;
for i3 = 0:(max_level_3 - 1)
    Z3 = Z3 + stat_weight_mode(i3, 3) * exp(-i3 * eps_3 / (AD.k * T_val));
end

rate_sum = 0;
k10_val = 0;
for i3 = 1:(max_level_3 - 1)
    x_i3 = stat_weight_mode(i3, 3) * exp(-i3 * eps_3 / (AD.k * T_val)) / Z3;
    if x_i3 < 1e-15
        continue
    end

    state_i = rxn.state_i;
    state_f = rxn.state_f;
    state_i(3) = i3;
    state_f(3) = i3 - 1;

    k_if = rates_fho_vv(T_val, state_i, state_f, rxn.state_k, rxn.state_kf, ...
        rxn.steric_vv, alpha, e_m, 'g0_max', 80.0, 'n_steps', 10000, ...
        'gamma', rxn.gamma);

    rate_sum = rate_sum + x_i3 * k_if;
    if i3 == 1
        k10_val = k_if;
    end
end

if k10_val > 0
    ptau_simpl = AD.k * T_val / k10_val / ATM;
else
    ptau_simpl = Inf;
end

c_vib = vibrational_heat_capacity_local(T_val, AD);
dim_eps = eps_3 / (AD.k * T_val);
rhs = (2.0 * pi * AD.k) / (AD.m * c_vib) * dim_eps^2 * rate_sum;

if rhs > 0 && isfinite(rhs)
    ptau_kin = AD.k * T_val / rhs / ATM;
else
    ptau_kin = Inf;
end
end

function c_vib = vibrational_heat_capacity_local(T, AD)
eps = AD.e1234 - AD.e1234(1);
beta = eps / (AD.k * T);
boltz = AD.stw(:).' .* exp(-beta(:).');
Z = sum(boltz);

S1 = sum(boltz .* beta(:).') / Z;
S2 = sum(boltz .* beta(:).'.^2) / Z;
c_vib = (AD.k / AD.m) * (S2 - S1^2);
end

function sw = stat_weight_mode(level, mode)
switch mode
    case 1
        sw = 1;
    case 2
        sw = level + 1;
    case 3
        sw = 0.5 * (level + 1) * (level + 2);
    case 4
        sw = 0.5 * (level + 1) * (level + 2);
    otherwise
        error('Invalid mode: %d', mode);
end
end

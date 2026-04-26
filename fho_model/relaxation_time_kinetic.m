function result = relaxation_time_kinetic(T, alpha, e_m, steric2, steric4, varargin)
% RELAXATION_TIME_KINETIC  VT & VV relaxation times using kinetic theory.
%
%   result = relaxation_time_kinetic(T, alpha, e_m, steric2, steric4)
%   result = relaxation_time_kinetic(..., 'vv_cases', {case1; case2; ...})
%
%   Computes p*tau [atm*s] using the kinetic theory formulation:
%
%     kT / (p*tau) = (2*pi*k_B) / (m * c_vibr) * (delta_eps/(kT))^2
%                    * sum_{i >= 1} x(i) * k_{i -> i-1}(T)
%
%   VT:  delta_eps = hc * omega_m  (one quantum of VT mode m = 2 or 4)
%        Sum over mode-m Boltzmann levels, rates from rates_fho_vt.
%
%   VV:  delta_eps = hc * omega_3  (one quantum of donor mode 3)
%        Sum over mode-3 Boltzmann levels, rates from rates_fho_vv.
%        (see cell 95 for derivation)
%
%   Inputs:
%     T       - temperature [K] (scalar or vector)
%     alpha   - common Morse range parameter [m^-1]
%     e_m     - common Morse well depth / k_B [K]
%     steric2 - steric factor for mode 2
%     steric4 - steric factor for mode 4
%
%   Optional name-value pairs:
%     'n_steps'       - integration points for rate calculation (default: 10000)
%     'pop_threshold' - minimum population fraction threshold (default: 1e-15)
%     'vv_cases'      - cell array of VV processes, each row:
%                        {si0, sf0, sk, skf, svv, label}
%                        or {si0, sf0, sk, skf, svv, gamma, label}
%                        si0, sf0: template states (mode 3 = 1->0)
%                        sk, skf:  partner molecule states
%                        svv: VV steric factor
%                        gamma: VV coupling mass-ratio parameter
%                        label: process label string
%
%   Output:
%     result  - struct with fields:
%       .ptau2      - p*tau for VT mode 2 [atm*s]
%       .ptau4      - p*tau for VT mode 4 [atm*s]
%       .ptau_total - total VT relaxation time [atm*s]
%       .ptau_vv    - (if vv_cases given) cell array of VV p*tau vectors
%       .ptau_vv_simpl - (if vv_cases given) simplified kT/k10 for comparison
%       .labels_vv  - (if vv_cases given) labels for each VV process
%
%   Requires: rates_fho_vt.m, rates_fho_vv.m, input_data.m

% --- parse optional arguments ---
p = inputParser;
addParameter(p, 'n_steps',       10000, @isnumeric);
addParameter(p, 'pop_threshold', 1e-15, @isnumeric);
addParameter(p, 'vv_cases',      {},    @iscell);
parse(p, varargin{:});
n_steps       = p.Results.n_steps;
pop_threshold = p.Results.pop_threshold;
vv_cases      = p.Results.vv_cases;

AD = input_data();

k_B     = AD.k;
mass    = AD.m;
levels  = AD.lch4;          % [9, 17, 9, 20]
omega   = AD.omega;         % [cm^-1]
hc      = AD.h * AD.c;      % h*c [J*cm]

% ── VT relaxation times ──
N = length(T);
ptau2     = zeros(1, N);
ptau4     = zeros(1, N);
ptau_tot  = zeros(1, N);

for iT = 1:N
    Ti = T(iT);
    c_vib = vibrational_heat_capacity(Ti, AD);

    ptau2(iT) = compute_ptau_vt_mode(Ti, 2, alpha, e_m, steric2, ...
                                     levels, omega, hc, k_B, mass, c_vib, ...
                                     n_steps, pop_threshold);

    ptau4(iT) = compute_ptau_vt_mode(Ti, 4, alpha, e_m, steric4, ...
                                     levels, omega, hc, k_B, mass, c_vib, ...
                                     n_steps, pop_threshold);

    alpha_pop    = relative_boltzmann_population_full(Ti, AD);
    ptau_tot(iT) = (1 + alpha_pop) / (1/ptau4(iT) + alpha_pop/ptau2(iT));
end

result.ptau2      = ptau2;
result.ptau4      = ptau4;
result.ptau_total = ptau_tot;

% ── VV relaxation times (kinetic theory) ──
if ~isempty(vv_cases)
    n_vv = size(vv_cases, 1);
    ptau_vv_kin   = cell(n_vv, 1);
    ptau_vv_simpl = cell(n_vv, 1);
    labels_vv     = cell(n_vv, 1);

    % delta_eps for VV = one quantum of donor mode 3
    eps_3       = hc * omega(3);            % [J]
    max_level_3 = levels(3);                % 9

    for ic = 1:n_vv
        si0      = vv_cases{ic, 1};
        sf0      = vv_cases{ic, 2};
        sk       = vv_cases{ic, 3};
        skf      = vv_cases{ic, 4};
        svv      = vv_cases{ic, 5};
        if size(vv_cases, 2) >= 7
            gamma_vv = vv_cases{ic, 6};
            labels_vv{ic} = vv_cases{ic, 7};
        else
            gamma_vv = 0.5;
            labels_vv{ic} = vv_cases{ic, 6};
        end

        ptau_vv_kin{ic}   = zeros(1, N);
        ptau_vv_simpl{ic} = zeros(1, N);

        for iT = 1:N
            Ti    = T(iT);
            c_vib = vibrational_heat_capacity(Ti, AD);

            [ptau_k, ptau_s] = compute_ptau_vv_mode(Ti, si0, sf0, sk, skf, ...
                                                     svv, alpha, e_m, gamma_vv, ...
                                                     eps_3, max_level_3, ...
                                                     k_B, mass, c_vib, ...
                                                     n_steps, pop_threshold);
            ptau_vv_kin{ic}(iT)   = ptau_k;
            ptau_vv_simpl{ic}(iT) = ptau_s;
        end
    end

    result.ptau_vv       = ptau_vv_kin;
    result.ptau_vv_simpl = ptau_vv_simpl;
    result.labels_vv     = labels_vv;
end

end


% =====================================================================
% VT: mode-specific p*tau using kinetic theory
% =====================================================================
function ptau = compute_ptau_vt_mode(T, mode, alpha, e_m, steric, ...
                                     levels, omega, hc, k_B, mass, c_vib, ...
                                     n_steps, pop_threshold)

if mode == 2
    mode_idx  = 2;
    max_level = levels(2);
elseif mode == 4
    mode_idx  = 4;
    max_level = levels(4);
else
    error('Only modes 2 and 4 supported, got %d', mode);
end

eps_1 = hc * omega(mode_idx);   % one quantum energy [J]

stat_w = zeros(1, max_level);
for i_m = 0:max_level-1
    stat_w(i_m+1) = stat_weight_mode(i_m, mode);
end
boltz = stat_w .* exp(-(0:max_level-1) * eps_1 / (k_B * T));
Z_m   = sum(boltz);

rate_sum = 0;
for i_m = 1:max_level-1
    x_im = stat_w(i_m+1) * exp(-i_m * eps_1 / (k_B * T)) / Z_m;
    if x_im < pop_threshold
        continue
    end
    k_if = rates_fho_vt(T, i_m, i_m - 1, mode_idx, steric, alpha, e_m, ...
                        'n_steps', n_steps);
    rate_sum = rate_sum + x_im * k_if;
end

if rate_sum <= 0
    ptau = 1e-20;
    return
end

dim_eps   = eps_1 / (k_B * T);
rhs       = (2 * pi * k_B) / (mass * c_vib) * dim_eps^2 * rate_sum;

if rhs <= 0 || ~isfinite(rhs)
    ptau = 1e-20;
    return
end

ptau = k_B * T / rhs / 101325;   % [atm*s]

end


% =====================================================================
% VV: kinetic theory p*tau for a single VV process
%   delta_eps = eps_3 = hc * omega_3  (one quantum of donor mode 3)
%   Sum runs over mode-3 Boltzmann levels
% =====================================================================
function [ptau_kin, ptau_simpl] = compute_ptau_vv_mode(T, si0, sf0, sk, skf, ...
                                                        steric_vv, ...
                                                        alpha_vv, e_m_vv, gamma_vv, ...
                                                        eps_3, max_level_3, ...
                                                        k_B, mass, c_vib, ...
                                                        n_steps, pop_threshold)

% Mode-3 partition function
Z3 = 0;
for i3 = 0:max_level_3-1
    Z3 = Z3 + stat_weight_mode(i3, 3) * exp(-i3 * eps_3 / (k_B * T));
end

rate_sum = 0;
k10_val  = 0;

for i3 = 1:max_level_3-1
    x_i3 = stat_weight_mode(i3, 3) * exp(-i3 * eps_3 / (k_B * T)) / Z3;
    if x_i3 < pop_threshold
        continue
    end

    si_lev    = si0;
    sf_lev    = sf0;
    si_lev(3) = i3;
    sf_lev(3) = i3 - 1;

    k_if = rates_fho_vv(T, si_lev, sf_lev, sk, skf, ...
                        steric_vv, alpha_vv, e_m_vv, ...
                        'n_steps', n_steps, 'gamma', gamma_vv);

    rate_sum = rate_sum + x_i3 * k_if;

    if i3 == 1
        k10_val = k_if;
    end
end

% Kinetic theory
dim_eps = eps_3 / (k_B * T);
rhs     = (2 * pi * k_B) / (mass * c_vib) * dim_eps^2 * rate_sum;

if rhs > 0 && isfinite(rhs)
    ptau_kin = k_B * T / rhs / 101325;
else
    ptau_kin = Inf;
end

% Simplified: kT / k10
if k10_val > 0
    ptau_simpl = k_B * T / k10_val / 101325;
else
    ptau_simpl = Inf;
end

end


% =====================================================================
% Vibrational heat capacity c_vibr(T) [J/(kg*K)]
% =====================================================================
function c_vib = vibrational_heat_capacity(T, AD)
% c_vibr = (k_B / m) * [<(eps/(kT))^2> - <eps/(kT)>^2]

eps   = AD.e1234 - AD.e1234(1);   % energy from ground state [J]
beta  = eps / (AD.k * T);
boltz = AD.stw(:).' .* exp(-beta(:).');
Z     = sum(boltz);

S1 = sum(boltz .* beta(:).') / Z;
S2 = sum(boltz .* beta(:).'.^2) / Z;

c_vib = (AD.k / AD.m) * (S2 - S1^2);

end


% =====================================================================
% Relative Boltzmann population of first excited mode 2 to mode 4 levels
% =====================================================================
function alpha_pop = relative_boltzmann_population_full(T, AD)
% alpha = (s_2(1) / s_4(1)) * exp(-(theta_2 - theta_4) / T)
% where theta_m = h*c*omega_m/k_B.

hc    = AD.h * AD.c;
omega = AD.omega;
k_B   = AD.k;

theta2 = hc * omega(2) / k_B;
theta4 = hc * omega(4) / k_B;
stat_ratio = stat_weight_mode(1, 2) / stat_weight_mode(1, 4);

alpha_pop = stat_ratio .* exp(-(theta2 - theta4) ./ T);

end


% =====================================================================
% Statistical weight for level n in mode m
% =====================================================================
function sw = stat_weight_mode(level, mode)
% Statistical weight for a given vibrational level and mode of CH4.

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

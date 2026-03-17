function result = relaxation_time_kinetic(T, alpha, e_m, steric2, steric4, varargin)
% RELAXATION_TIME_KINETIC  VT relaxation times using kinetic theory formulation.
%
%   result = relaxation_time_kinetic(T, alpha, e_m, steric2, steric4)
%   result = relaxation_time_kinetic(T, alpha, e_m, steric2, steric4, 'n_steps', 10000)
%
%   Computes p*tau [atm*s] for VT2, VT4 and the total VT relaxation time 
%   using the kinetic theory based formulation:
%
%     kT / (p*tau_m) = (2*pi*k_B) / (m * c_vibr) * (delta_eps/(kT))^2 
%                      * sum_{i_m >= 1} x_m(i_m) * k_{i_m -> i_m-1}(T)
%
%   where:
%     c_vibr(T)        - full vibrational heat capacity [J/(kg*K)]
%     x_m(i_m)         - molar fraction for level i_m in mode m
%     k_{i -> i-1}(T)  - FHO VT rate coefficient for i_m -> i_m-1
%     delta_eps        - one quantum energy of mode m
%
%   Inputs:
%     T       - temperature [K] (scalar or vector)
%     alpha   - Morse potential range parameter [m^-1]
%     e_m     - Morse well depth / k_B [K]
%     steric2 - steric factor for mode 2
%     steric4 - steric factor for mode 4
%
%   Optional name-value pairs:
%     'n_steps'       - integration points for rate calculation (default: 10000)
%     'pop_threshold' - minimum population fraction threshold (default: 1e-15)
%
%   Output:
%     result  - struct with fields:
%       .ptau2   - p*tau for mode 2 [atm*s] (same size as T)
%       .ptau4   - p*tau for mode 4 [atm*s] (same size as T)
%       .ptau_total - total VT relaxation time [atm*s] (same size as T)
%
%   Requires: rates_fho_vt.m, probabilties_fho_vt.m, input_data.m, vdf.m

% --- parse optional arguments ---
p = inputParser;
addParameter(p, 'n_steps',       10000, @isnumeric);
addParameter(p, 'pop_threshold', 1e-15, @isnumeric);
parse(p, varargin{:});
n_steps       = p.Results.n_steps;
pop_threshold = p.Results.pop_threshold;

AD = input_data();

k_B     = AD.k;
mass    = AD.m;
levels  = AD.lch4;          % [9, 17, 9, 20]
d_modes = AD.d;             % [1, 2, 3, 3]
omega   = AD.omega;         % [cm^-1]
hc      = AD.h * AD.c;      % h*c [J*cm]

% pre-allocate outputs
N = length(T);
ptau2     = zeros(1, N);
ptau4     = zeros(1, N);
ptau_tot  = zeros(1, N);

for iT = 1:N
    Ti = T(iT);

    % --- vibrational heat capacity ---
    c_vib = vibrational_heat_capacity(Ti, AD);

    % --- mode 2: p*tau_2 ---
    ptau2(iT) = compute_ptau_mode(Ti, 2, alpha, e_m, steric2, ...
                                  levels, omega, hc, k_B, mass, c_vib, ...
                                  n_steps, pop_threshold);

    % --- mode 4: p*tau_4 ---
    ptau4(iT) = compute_ptau_mode(Ti, 4, alpha, e_m, steric4, ...
                                  levels, omega, hc, k_B, mass, c_vib, ...
                                  n_steps, pop_threshold);

    % --- total VT relaxation time ---
    alpha_pop   = relative_boltzmann_population_full(Ti, AD);
    ptau_tot(iT) = (1 + alpha_pop) / (1/ptau4(iT) + alpha_pop/ptau2(iT));
end

result.ptau2      = ptau2;
result.ptau4      = ptau4;
result.ptau_total = ptau_tot;

end


% =====================================================================
% Mode-specific p*tau using kinetic theory
% =====================================================================
function ptau = compute_ptau_mode(T, mode, alpha, e_m, steric, ...
                                  levels, omega, hc, k_B, mass, c_vib, ...
                                  n_steps, pop_threshold)
% Compute p*tau [atm*s] for a single mode.

if mode == 2
    mode_idx  = 2;      % MATLAB 1-based
    max_level = levels(2);
elseif mode == 4
    mode_idx  = 4;
    max_level = levels(4);
else
    error('Only modes 2 and 4 supported, got %d', mode);
end

% one quantum energy [J]
eps_1 = hc * omega(mode_idx);

% --- single-mode partition function ---
stat_w = zeros(1, max_level);
for i_m = 0:max_level-1
    stat_w(i_m+1) = stat_weight_mode(i_m, mode);
end
boltz = stat_w .* exp(-(0:max_level-1) * eps_1 / (k_B * T));
Z_m   = sum(boltz);

% --- weighted rate sum ---
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

dim_eps = eps_1 / (k_B * T);

% kT/(p*tau) = (2*pi*k_B)/(m*c_vib) * (eps/(kT))^2 * rate_sum  [m^3/s]
rhs = (2 * pi * k_B) / (mass * c_vib) * dim_eps^2 * rate_sum;

if rhs <= 0 || ~isfinite(rhs)
    ptau = 1e-20;
    return
end

ptau_Pa_s = k_B * T / rhs;       % [Pa*s]
ptau      = ptau_Pa_s / 101325;   % [atm*s]

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
% Relative Boltzmann population (energy-weighted) of mode 2 to mode 4
% =====================================================================
function alpha_pop = relative_boltzmann_population_full(T, AD)
% alpha = <E_2> / <E_4> where averages are over single-mode distributions

hc    = AD.h * AD.c;
omega = AD.omega;
k_B   = AD.k;

levels2 = AD.lch4(2);
levels4 = AD.lch4(4);

% mode 2
sw2 = zeros(1, levels2);
e2  = zeros(1, levels2);
for n = 0:levels2-1
    sw2(n+1) = stat_weight_mode(n, 2);
    e2(n+1)  = hc * omega(2) * n;
end
b2  = sw2 .* exp(-e2 / (k_B * T));
Z2  = sum(b2);
E2  = sum(e2 .* b2) / Z2;

% mode 4
sw4 = zeros(1, levels4);
e4  = zeros(1, levels4);
for n = 0:levels4-1
    sw4(n+1) = stat_weight_mode(n, 4);
    e4(n+1)  = hc * omega(4) * n;
end
b4  = sw4 .* exp(-e4 / (k_B * T));
Z4  = sum(b4);
E4  = sum(e4 .* b4) / Z4;

alpha_pop = E2 / E4;

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

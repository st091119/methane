function ptau = relaxation_time_kinetic_vv(T, si0, sf0, sk, skf, ...
                                            steric_vv, alpha, e_m, varargin)
% RELAXATION_TIME_KINETIC_VV  VV relaxation time using kinetic theory formulation.
% ! Simplified evaluation
%
%   ptau = relaxation_time_kinetic_vv(T, si0, sf0, sk, skf,
%                                     steric_vv, alpha, e_m)
%
%   Computes p*tau [atm*s] for a VV process using the kinetic theory
%   formulation (same structure as VT kinetic):
%
%     kT / (p*tau) = (2*pi*k_B) / (m * c_vibr) * (eps_donor/(kT))^2
%                    * sum_{i3 >= 1} x_3(i3) * k_{i3 -> i3-1}(T)
%
%   The donor mode is always mode 3 (nu_3). The summation goes over
%   mode-3 vibrational levels with Boltzmann weights.
%
%   Inputs:
%     T          - temperature [K] (scalar or vector)
%     si0        - template initial state of mol 1 [1x4], mode 3 = 1
%     sf0        - template final state of mol 1   [1x4], mode 3 = 0
%     sk         - initial state of mol 2 [1x4]
%     skf        - final state of mol 2   [1x4]
%     steric_vv  - steric factor S_VV
%     alpha      - Morse potential range parameter [m^-1]
%     e_m        - Morse well depth / k_B [K]
%
%   Optional name-value pairs:
%     'n_steps'       - integration points for rate (default: 10000)
%     'pop_threshold' - minimum population fraction (default: 1e-15)
%     'gamma'         - VV coupling mass-ratio parameter (default: 0.5)
%
%   Output:
%     ptau  - p*tau [atm*s], same size as T
%
%   Requires: rates_fho_vv.m, input_data.m

% --- parse optional arguments ---
p = inputParser;
addParameter(p, 'n_steps',       10000, @isnumeric);
addParameter(p, 'pop_threshold', 1e-15, @isnumeric);
addParameter(p, 'gamma',         0.5,   @isnumeric);
parse(p, varargin{:});
n_steps       = p.Results.n_steps;
pop_threshold = p.Results.pop_threshold;
gamma         = p.Results.gamma;

AD = input_data();

k_B     = AD.k;
mass    = AD.m;
levels  = AD.lch4;       % [9, 17, 9, 20]
omega   = AD.omega;       % [cm^-1]
hc      = AD.h * AD.c;    % h*c [J*cm]

% Donor mode = mode 3 (MATLAB index 3)
donor_mode  = 3;
max_level_3 = levels(donor_mode);   % 9
eps_3       = hc * omega(donor_mode);  % one quantum [J]

N    = length(T);
ptau = zeros(1, N);

for iT = 1:N
    Ti = T(iT);

    % --- mode-3 partition function ---
    Z3 = 0;
    for i3 = 0:max_level_3-1
        sw = stat_weight_mode(i3, 3);
        Z3 = Z3 + sw * exp(-i3 * eps_3 / (k_B * Ti));
    end

    % --- weighted rate sum ---
    rate_sum = 0;
    for i3 = 1:max_level_3-1
        x_i3 = stat_weight_mode(i3, 3) * exp(-i3 * eps_3 / (k_B * Ti)) / Z3;
        if x_i3 < pop_threshold
            continue
        end

        si_lev    = si0;
        sf_lev    = sf0;
        si_lev(3) = i3;        % mode 3 of initial state
        sf_lev(3) = i3 - 1;    % mode 3 of final state

        k_if = rates_fho_vv(Ti, si_lev, sf_lev, sk, skf, ...
                            steric_vv, alpha, e_m, ...
                            'n_steps', n_steps, 'gamma', gamma);
        rate_sum = rate_sum + x_i3 * k_if;
    end

    % --- kinetic theory ---
    c_vib   = vibrational_heat_capacity(Ti, AD);
    dim_eps = eps_3 / (k_B * Ti);
    rhs     = (2 * pi * k_B) / (mass * c_vib) * dim_eps^2 * rate_sum;

    if rhs > 0 && isfinite(rhs)
        ptau_Pa  = k_B * Ti / rhs;
        ptau(iT) = ptau_Pa / 101325;
    else
        ptau(iT) = Inf;
    end
end

end


% =====================================================================
% Vibrational heat capacity c_vibr(T) [J/(kg*K)]
% =====================================================================
function c_vib = vibrational_heat_capacity(T, AD)

eps   = AD.e1234 - AD.e1234(1);
beta  = eps / (AD.k * T);
boltz = AD.stw(:).' .* exp(-beta(:).');
Z     = sum(boltz);

S1 = sum(boltz .* beta(:).') / Z;
S2 = sum(boltz .* beta(:).'.^2) / Z;

c_vib = (AD.k / AD.m) * (S2 - S1^2);

end


% =====================================================================
% Statistical weight for level n in mode m
% =====================================================================
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

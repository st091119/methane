function k_VV = rates_fho_vv(T, state_i, state_f, state_k, state_kf, ...
                              steric_vv, alpha, e_m, varargin)
% RATES_FHO_VV  Thermally-averaged VV rate coefficient k_VV(T) [m^3/s].
%
%   k_VV = rates_fho_vv(T, state_i, state_f, state_k, state_kf,
%                        steric_vv, alpha, e_m)
%   k_VV = rates_fho_vv(..., 'g0_max', 80, 'n_steps', 10000)
%
%   Computes the thermally-averaged VV rate coefficient by integrating
%   the product sigma_coll * P_VV over a Maxwellian velocity distribution:
%
%     k = 4 * sqrt(2 k T / (pi m_r)) * int_0^inf exp(-g0^2) g0^3 sigma P_VV dg0
%
%   where g0 = g / sqrt(2 k T / m_r).
%
%   Inputs:
%     T          - temperature [K] (scalar)
%     state_i    - [i1 i2 i3 i4] initial state, molecule 1
%     state_f    - [f1 f2 f3 f4] final state,   molecule 1
%     state_k    - [k1 k2 k3 k4] initial state, molecule 2
%     state_kf   - [k1 k2 k3 k4] final state,   molecule 2
%     steric_vv  - steric factor S_VV (enters rho)
%     alpha      - Morse potential range parameter [m^-1]
%     e_m        - Morse well depth / k_B [K] (kept for parameter-set compatibility)
%
%   Optional name-value pairs:
%     'g0_max'  - upper integration limit (default: 80)
%     'n_steps' - number of grid points   (default: 10000)
%     'gamma'   - VV coupling mass-ratio parameter (default: 0.5)
%
%   Output:
%     k_VV  - rate coefficient [m^3/s]
%
%   Requires: probabilties_fho_vv.m, input_data.m

% --- parse optional arguments ---
p = inputParser;
addParameter(p, 'g0_max',  80,    @isnumeric);
addParameter(p, 'n_steps', 10000, @isnumeric);
addParameter(p, 'gamma',   0.5,   @isnumeric);
parse(p, varargin{:});
g0_max  = p.Results.g0_max;
n_steps = p.Results.n_steps;
gamma   = p.Results.gamma;

AD = input_data();

% --- collision parameters ---
kk   = AD.k;            % Boltzmann constant [J/K]
mass = AD.m;             % CH4 mass [kg]
m_r  = mass / 2;         % reduced mass [kg]

% --- VSS parameters for CH4-CH4 ---
omega_d   = 0.185521;
omega_eta = 0.174955;
c_d       = 1.492958e-18;
c_eta     = 1.074084e-18;

% --- dimensionless velocity grid ---
g0 = linspace(1e-5, g0_max, n_steps);
g  = g0 * sqrt(2 * kk * T / m_r);

% --- VV transition probability ---
P_VV = probabilties_fho_vv(g, state_i, state_f, state_k, state_kf, ...
                            steric_vv, alpha, e_m, gamma);

% --- VSS collision cross section sigma_coll(g) [m^2] ---
x        = m_r * g.^2 / (2 * kk);
sigma_d  = c_d   .* x .^ (-omega_d);
sigma_e  = c_eta .* x .^ (-omega_eta);
sigma    = 0.5 .* sigma_d .* (2 * sigma_d + sigma_e) ./ (2 * sigma_d - sigma_e);

% --- integrand ---
integrand = exp(-g0.^2) .* g0.^3 .* sigma .* P_VV;

% --- trapezoidal integration with NaN/Inf protection ---
integrand(~isfinite(integrand)) = 0;
integral_val = trapz(g0, integrand);

% --- prefactor ---
pref = 4 * sqrt(2 * kk * T / (pi * m_r));

k_VV = pref * integral_val;

end

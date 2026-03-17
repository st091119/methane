function k_VT = rates_fho_vt(T, qi, qf, mode, steric, alpha, e_m, varargin)
% RATES_FHO_VT  Thermally-averaged VT rate coefficient k_VT(T) [m^3/s].
%
%   k_VT = rates_fho_vt(T, qi, qf, mode, steric, alpha, e_m)
%   k_VT = rates_fho_vt(T, qi, qf, mode, steric, alpha, e_m, 'g0_max', 80, 'n_steps', 10000)
%
%   Computes the thermally-averaged VT rate coefficient by integrating
%   the product sigma_coll * P_VT over a Maxwellian velocity distribution:
%
%     k = 4 * sqrt(2 k T / (pi m_r)) * int_0^inf exp(-g0^2) g0^3 sigma_coll P_VT dg0
%
%   where g0 = g / sqrt(2 k T / m_r) is the dimensionless relative velocity.
%
%   Inputs:
%     T       - temperature [K] (scalar)
%     qi      - initial quantum number of the active mode
%     qf      - final quantum number of the active mode
%     mode    - active mode index: 1, 2, 3, or 4
%     steric  - steric factor S_VT
%     alpha   - Morse potential range parameter [m^-1]
%     e_m     - Morse well depth / k_B [K]
%
%   Optional name-value pairs:
%     'g0_max'  - upper integration limit in g0 (default: 80)
%     'n_steps' - number of integration points (default: 10000)
%
%   Output:
%     k_VT    - rate coefficient [m^3/s]
%
%   Requires: probabilties_fho_vt.m, input_data.m

% --- parse optional arguments ---
p = inputParser;
addParameter(p, 'g0_max',  80,    @isnumeric);
addParameter(p, 'n_steps', 10000, @isnumeric);
parse(p, varargin{:});
g0_max  = p.Results.g0_max;
n_steps = p.Results.n_steps;

AD = input_data();

% --- collision parameters ---
k    = AD.k;             % Boltzmann constant [J/K]
mass = AD.m;             % CH4 mass [kg]
m_r  = mass / 2;         % reduced mass [kg]

% --- VSS parameters for CH4-CH4 ---
omega_d   = 0.185521;
omega_eta = 0.174955;
c_d       = 1.492958e-18;
c_eta     = 1.074084e-18;

% --- dimensionless velocity grid ---
g0 = linspace(1e-5, g0_max, n_steps);
g  = g0 * sqrt(2 * k * T / m_r);

% --- VT transition probability ---
P_VT = probabilties_fho_vt(g, qi, qf, mode, steric, alpha, e_m);

% --- VSS collision cross section sigma_coll(g) [m^2] ---
x        = m_r * g.^2 / (2 * k);
sigma_d  = c_d   .* x .^ (-omega_d);
sigma_e  = c_eta .* x .^ (-omega_eta);
sigma    = 0.5 .* sigma_d .* (2 * sigma_d + sigma_e) ./ (2 * sigma_d - sigma_e);

% --- integrand ---
integrand = exp(-g0.^2) .* g0.^3 .* sigma .* P_VT;

% --- trapezoidal integration with NaN/Inf protection ---
integrand(~isfinite(integrand)) = 0;
integral_val = trapz(g0, integrand);

% --- prefactor ---
pref = 4 * sqrt(2 * k * T / (pi * m_r));

k_VT = pref * integral_val;

end

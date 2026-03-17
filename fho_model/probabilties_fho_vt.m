function P = probabilties_fho_vt(g, qi, qf, mode, steric, alpha, e_m)
% PROBABILTIES_FHO_VT  Original (factorial) FHO VT transition probability for only CH4-CH4 collisions.
%
%   P = probabilties_fho_vt(g, qi, qf, mode, steric, alpha, e_m)
%
%   Computes the FHO VT transition probability P_VT(g) for a single
%   active mode of CH4 using the original factorial formulation:
%
%     P = qi! * qf! * eps^(qi+qf) * exp(-eps) 
%         * [ sum_{r=0}^{min(qi,qf)} (-1)^r / (r!(qi-r)!(qf-r)! eps^r) ]^2
%
%   Inputs:
%     g       - relative collision velocities [m/s], row or column vector
%     qi      - initial quantum number of the active mode (integer >= 0)
%     qf      - final quantum number of the active mode (integer >= 0)
%     mode    - active mode index: 1, 2, 3, or 4
%     steric  - steric factor S_VT (dimensionless)
%     alpha   - Morse potential range parameter [m^-1]
%     e_m     - Morse well depth / k_B [K]
%
%   Output:
%     P       - VT transition probability (same size as g), clamped to [0, 1]
%
%   Requires: input_data.m (for oscillator reduced masses and mode frequencies)

g = g(:).';  % ensure row vector

AD = input_data();

% --- physical constants ---
hbar = AD.h / (2*pi);       % reduced Planck constant [J*s]
k    = AD.k;                % Boltzmann constant [J/K]
c_SI = 299792458;            % speed of light [m/s]

% --- oscillator reduced masses [kg] 
mu_osc = [1.6734911e-27, 1.6734911e-27, 1.8305781e-27, 1.9562809e-27];

% --- mode frequencies [m^-1] ---
omega_m_e = [302550, 158270, 315680, 136740];

% --- collision parameters ---
mass  = AD.m;                % CH4 mass [kg]
m_r   = mass / 2;            % reduced mass [kg]
gamma = 0.5;                 % mass ratio parameter

% --- total quanta exchanged ---
s = abs(qi - qf);
if s == 0
    P = zeros(size(g));
    return
end

% --- energy gap |delta_E| [J] ---
delta_e_abs = AD.h * c_SI * omega_m_e(mode) * s;

% --- adiabaticity parameter epsilon_VT ---
eps = epsilon_vt_calc(g, delta_e_abs, s, mu_osc(mode), steric, ...
                      m_r, gamma, alpha, e_m, hbar, k);

% --- Laguerre-type sum ---
n_min = min(qi, qf);

% pre-compute coefficients
coeffs = zeros(1, n_min + 1);
for r = 0:n_min
    coeffs(r+1) = ((-1)^r) / (factorial(r) * factorial(qi - r) * factorial(qf - r));
end

% inverse powers of epsilon
eps_safe = max(eps, 1e-300);   % avoid division by zero
inv_eps_pow = ones(n_min + 1, length(g));
for r = 1:n_min
    inv_eps_pow(r+1, :) = inv_eps_pow(r, :) ./ eps_safe;
end

% sum over r
laguerre_sum = zeros(1, length(g));
for r = 0:n_min
    laguerre_sum = laguerre_sum + coeffs(r+1) * inv_eps_pow(r+1, :);
end

% --- full probability ---
prefactor = factorial(qi) * factorial(qf);
P = prefactor .* (eps .^ (qi + qf)) .* exp(-eps) .* laguerre_sum.^2;

% zero out where epsilon was zero
P(eps <= 0) = 0;

% clamp to [0, 1]
P = min(max(P, 0), 1);

end


% =====================================================================
% Adiabaticity parameter epsilon_VT
% =====================================================================
function eps = epsilon_vt_calc(g, delta_e_abs, s, mu_osc_m, steric, ...
                               m_r, gamma, alpha, e_m, hbar, k)
% EPSILON_VT_CALC  Single-quantum VT adiabaticity parameter.
%
%   eps = (4 pi omega gamma^2 m_r^2) / (hbar alpha^2 mu_osc) 
%         * cosh^2(a) / sinh^2(b) * S_steric
%
%   where:
%     omega = |delta_E| / (hbar * s)
%     a = (1 + phi) * pi * omega / (alpha * v_bar)
%     b = 2 * pi * omega / (alpha * v_bar)

omega = delta_e_abs / (hbar * s);

% symmetrised velocity
v2   = sqrt(g.^2 + delta_e_abs / m_r);
vbar = 0.5 * (v2 + g);

% phase angle correction
phi = (2/pi) * atan(sqrt(2 * (e_m * k) ./ (m_r * vbar.^2)));

a = (1 + phi) .* pi .* omega ./ (alpha .* vbar);
b = 2 * pi * omega ./ (alpha .* vbar);

% safe cosh^2 / sinh^2
ratio = safe_cosh2_over_sinh2(a, b);

pref = (4 / hbar) * (pi * omega * gamma^2 * m_r^2) / (alpha^2 * mu_osc_m);
eps  = pref .* ratio .* steric;

end


% =====================================================================
% Numerically-safe cosh^2(a) / sinh^2(b)
% =====================================================================
function result = safe_cosh2_over_sinh2(a, b)
% SAFE_COSH2_OVER_SINH2  Compute cosh^2(a) / sinh^2(b) with overflow protection.

a_abs = abs(a);
b_abs = abs(b);
result = zeros(size(a));

% both large
both_large = (a_abs > 50) & (b_abs > 50);
result(both_large) = exp(2 * (a_abs(both_large) - b_abs(both_large)));

% only a large
a_large = (a_abs > 50) & (b_abs <= 50);
if any(a_large)
    c2a = 0.25 * exp(2 * a_abs(a_large));
    s2b = sinh(b(a_large)).^2;
    s2b = max(s2b, 1e-100);
    result(a_large) = c2a ./ s2b;
end

% only b large
b_large = (a_abs <= 50) & (b_abs > 50);
if any(b_large)
    c2a = cosh(a(b_large)).^2;
    s2b = 0.25 * exp(2 * b_abs(b_large));
    result(b_large) = c2a ./ s2b;
end

% both moderate
both_small = (a_abs <= 50) & (b_abs <= 50);
if any(both_small)
    s2b = sinh(b(both_small)).^2;
    s2b = max(s2b, 1e-100);
    result(both_small) = cosh(a(both_small)).^2 ./ s2b;
end

end

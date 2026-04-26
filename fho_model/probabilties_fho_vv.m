function P = probabilties_fho_vv(g, state_i, state_f, state_k, state_kf, ...
                                  steric_vv, alpha, ~, gamma)
% PROBABILTIES_FHO_VV  Simplified FHO VV transition probability for CH4-CH4.
%
%   P = probabilties_fho_vv(g, state_i, state_f, state_k, state_kf,
%                           steric_vv, alpha, e_m, gamma)
%
%   This is the simplified VV model used by the Python FHO implementation.
%   The old Zelechow epsilon/VT-steric contribution is not used; only the VV
%   steric factor enters through the oscillator-oscillator coupling rho.
%   gamma defaults to the CH4 value 0.5; selected fitted VV channels can
%   override it through rates_fho_vv(..., 'gamma', value).

if nargin < 9
    gamma = 0.5;
end

g = g(:).';
state_i = state_i(:).';
state_f = state_f(:).';
state_k = state_k(:).';
state_kf = state_kf(:).';

AD = input_data();

hbar = AD.h / (2 * pi);
c_SI = 299792458;

mu_osc = [1.6734911e-27, 1.6734911e-27, 1.8305781e-27, 1.9562809e-27];
omega_m_e = [302550, 158270, 315680, 136740];

mass = AD.m;
m_r = mass / 2;

s1 = total_quanta(state_i, state_f);
s2 = total_quanta(state_k, state_kf);

if s1 == 0 && s2 == 0
    P = zeros(size(g));
    return
end

s = max(max(s1, s2), 1);

n1 = 1.0;
if s1 > 0
    n1 = n_s_factor(state_i, state_f, s1);
end

n2 = 1.0;
if s2 > 0
    n2 = n_s_factor(state_k, state_kf, s2);
end

n_s = n1 * n2;

delta1 = delta_e_state(state_i, state_f, AD.h, c_SI, omega_m_e);
delta2 = delta_e_state(state_k, state_kf, AD.h, c_SI, omega_m_e);
delta_tot_abs = abs(delta1 + delta2);

v2 = sqrt(g.^2 + delta_tot_abs / m_r);
vbar = 0.5 * (v2 + g);

if s1 > 0 && abs(delta1) > 0
    omega1 = abs(delta1) / (hbar * s1);
else
    omega1 = 0;
end

if s2 > 0 && abs(delta2) > 0
    omega2 = abs(delta2) / (hbar * s2);
else
    omega2 = 0;
end

idx1 = active_mode_index(state_i, state_f);
idx2 = active_mode_index(state_k, state_kf);
mu1 = mu_osc(idx1);
mu2 = mu_osc(idx2);

if omega1 == 0 && omega2 == 0
    P = zeros(size(g));
    return
elseif omega1 ~= 0 && omega2 ~= 0
    mu_eff = sqrt(mu1 * mu2);
    rho = 2 * (m_r / mu_eff) * gamma^2 * alpha .* vbar ...
        ./ sqrt(omega1 * omega2) * steric_vv;
elseif omega1 ~= 0
    rho = 2 * (m_r / mu1) * gamma^2 * alpha .* vbar ...
        ./ omega1 * steric_vv;
else
    rho = 2 * (m_r / mu2) * gamma^2 * alpha .* vbar ...
        ./ omega2 * steric_vv;
end

xi = (pi^2 ./ (4 * alpha .* vbar)) .* abs(omega1 - omega2);
rho_xi = rho .* safe_x_over_sinh(xi);
lambda_eff = (rho_xi.^2) / 4;

prefactor = n_s^s / (factorial(s)^2);
exponent = -2 * n_s .* lambda_eff / (s + 1);
P = prefactor .* (lambda_eff.^s) .* exp(exponent);
P(~isfinite(P)) = 0;
P = min(max(P, 0), 1);

end

function s = total_quanta(state_i, state_f)
s = sum(abs(state_i - state_f));
end

function n_s = n_s_factor(state_i, state_f, s)
if s == 0
    n_s = 1.0;
    return
end

prod_val = 1.0;
for m = 1:numel(state_i)
    imax = max(state_i(m), state_f(m));
    imin = min(state_i(m), state_f(m));
    prod_val = prod_val * factorial(imax) / factorial(imin);
end
n_s = prod_val^(1 / s);
end

function idx = active_mode_index(state_i, state_f)
[~, idx] = max(abs(state_f - state_i));
end

function de = delta_e_state(state_i, state_f, h, c_SI, omega_m_e)
de = h * c_SI * sum(omega_m_e .* (state_f - state_i));
end

function result = safe_x_over_sinh(x)
result = ones(size(x));

small = abs(x) < 1e-5;
large = abs(x) > 50;
moderate = ~(small | large);

if any(large)
    result(large) = 2 .* x(large) .* exp(-x(large));
end

if any(moderate)
    result(moderate) = x(moderate) ./ sinh(x(moderate));
end

result(~isfinite(result)) = 0;
end

function P = probabilties_fho_vv(g, state_i, state_f, state_k, state_kf, ...
                                  steric_vt, steric_vv, alpha, e_m)
% PROBABILTIES_FHO_VV  Original (Zelechow) FHO VV transition probability for CH4-CH4.
%
%   P = probabilties_fho_vv(g, state_i, state_f, state_k, state_kf,
%                           steric_vt, steric_vv, alpha, e_m)
%
%   Molecule 1 (donor):    state_i -> state_f
%   Molecule 2 (acceptor): state_k -> state_kf
%
%   Handles both intramolecular VV (state_k == state_kf, both modes on mol 1)
%   and intermolecular VV (both molecules change).
%
%   Inputs:
%     g          - relative collision velocities [m/s], vector
%     state_i    - [i1 i2 i3 i4] initial state, molecule 1
%     state_f    - [f1 f2 f3 f4] final state,   molecule 1
%     state_k    - [k1 k2 k3 k4] initial state, molecule 2
%     state_kf   - [k1 k2 k3 k4] final state,   molecule 2
%     steric_vt  - steric factor S_VT (enters epsilon)
%     steric_vv  - steric factor S_VV (enters rho)
%     alpha      - Morse potential range parameter [m^-1]
%     e_m        - Morse well depth / k_B [K]
%
%   Output:
%     P          - VV transition probability (same size as g), clamped to [0, 1]
%
%   Requires: input_data.m

g = g(:).';  % ensure row vector

AD = input_data();

% --- physical constants ---
hbar = AD.h / (2*pi);       % reduced Planck constant [J*s]
k    = AD.k;                % Boltzmann constant [J/K]
c_SI = 299792458;            % speed of light [m/s]

% --- oscillator reduced masses [kg] ---
mu_osc = [1.6734911e-27, 1.6734911e-27, 1.8305781e-27, 1.9562809e-27];

% --- mode frequencies [m^-1] ---
omega_m_e = [302550, 158270, 315680, 136740];

% --- collision parameters ---
mass  = AD.m;
m_r   = mass / 2;
gamma = 0.5;

% ── Check for intramolecular VV (molecule 2 unchanged) ──
mol2_unchanged = all(state_k == state_kf);

if mol2_unchanged
    % Intramolecular: donor and acceptor on the same molecule
    diff1 = state_f - state_i;   % per-mode change on mol 1

    losing_modes  = find(diff1 < 0);   % modes losing quanta
    gaining_modes = find(diff1 > 0);   % modes gaining quanta

    if isempty(losing_modes) || isempty(gaining_modes)
        P = zeros(size(g));
        return
    end

    % Donor = mode losing the most quanta
    [~, idx_max] = max(state_i(losing_modes) - state_f(losing_modes));
    donor_idx = losing_modes(idx_max);
    qi_donor  = state_i(donor_idx);
    qf_donor  = state_f(donor_idx);

    % Acceptor = gaining modes on mol 1
    qi_acceptor = sum(state_i(gaining_modes));
    qf_acceptor = sum(state_f(gaining_modes));

    % Energy contributions
    delta_donor_val    = AD.h * c_SI * sum(omega_m_e(losing_modes)  .* diff1(losing_modes));
    delta_acceptor_val = AD.h * c_SI * sum(omega_m_e(gaining_modes) .* diff1(gaining_modes));

    s_donor    = sum(state_i(losing_modes) - state_f(losing_modes));
    s_acceptor = sum(state_f(gaining_modes) - state_i(gaining_modes));

    omega_donor    = abs(delta_donor_val) / (hbar * s_donor);
    omega_acceptor = abs(delta_acceptor_val) / (hbar * s_acceptor);

    delta_tot_abs = abs(delta_donor_val + delta_acceptor_val);

    % Acceptor mode index (mode gaining the most quanta)
    [~, idx_acc_max] = max(state_f(gaining_modes) - state_i(gaining_modes));
    acceptor_idx = gaining_modes(idx_acc_max);

else
    % Standard intermolecular VV
    donor_diffs = abs(state_f - state_i);
    if max(donor_diffs) == 0
        P = zeros(size(g));
        return
    end
    [~, donor_idx] = max(donor_diffs);
    qi_donor = state_i(donor_idx);
    qf_donor = state_f(donor_idx);

    acceptor_diffs = abs(state_kf - state_k);
    acceptor_modes = find(acceptor_diffs > 0);
    if isempty(acceptor_modes)
        P = zeros(size(g));
        return
    end
    qi_acceptor = sum(state_k(acceptor_modes));
    qf_acceptor = sum(state_kf(acceptor_modes));

    % Energy gaps
    delta_donor    = AD.h * c_SI * sum(omega_m_e .* (state_f - state_i));
    delta_acceptor = AD.h * c_SI * sum(omega_m_e .* (state_kf - state_k));
    delta_tot_abs  = abs(delta_donor + delta_acceptor);

    s_donor    = sum(donor_diffs);
    s_acceptor = sum(acceptor_diffs);

    if abs(delta_donor) > 0
        omega_donor = abs(delta_donor) / (hbar * s_donor);
    else
        omega_donor = 0;
    end
    if abs(delta_acceptor) > 0
        omega_acceptor = abs(delta_acceptor) / (hbar * s_acceptor);
    else
        omega_acceptor = 0;
    end

    [~, acceptor_idx] = max(acceptor_diffs);
end

% ── Common path ──
i12 = qi_donor + qi_acceptor;
f12 = qf_donor + qf_acceptor;
n   = min(i12 + 1, f12 + 1);

if n == 0 || omega_donor == 0 || omega_acceptor == 0
    P = zeros(size(g));
    return
end

% ── Zelechow matrices ──
C_i12 = zelechow_C_matrix(i12);
C_f12 = zelechow_C_matrix(f12);

% ── Symmetrised velocity from total defect ──
v2   = sqrt(g.^2 + delta_tot_abs / m_r);
vbar = 0.5 * (v2 + g);

% ── Phase angle correction ──
phi = (2/pi) * atan(sqrt(2 * (e_m * k) ./ (m_r * vbar.^2)));

% ── epsilon parameter (donor mode, omega from donor, vbar from total defect) ──
a_arg = (1 + phi) .* pi .* omega_donor ./ (alpha .* vbar);
b_arg = 2 * pi * omega_donor ./ (alpha .* vbar);
ratio = safe_cosh2_over_sinh2(a_arg, b_arg);

pref_eps = (4 / hbar) * (pi * omega_donor * gamma^2 * m_r^2) / ...
           (alpha^2 * mu_osc(donor_idx));
eps = pref_eps .* ratio .* steric_vt;

% ── rho parameter (VV coupling) ──
mu_eff = sqrt(mu_osc(donor_idx) * mu_osc(acceptor_idx));
rho = 2 * (m_r / mu_eff) * gamma^2 * alpha .* vbar ...
      ./ sqrt(omega_donor * omega_acceptor) * steric_vv;

% ── Complex Zelechow sum (vectorised over g) ──
eps_safe = max(eps, 1e-300);
total_re = zeros(size(g));
total_im = zeros(size(g));

for j = 1:n
    a_val = i12 - j + 1;
    b_val = f12 - j + 1;

    c_i = C_i12(j, qi_acceptor + 1);    % 1-based indexing
    c_f = C_f12(j, qf_acceptor + 1);

    sgn       = (-1)^(i12 - j + 1);
    eps_pow   = eps_safe .^ ((a_val + b_val) / 2);
    exp_eps   = exp(-eps / 2);
    sqrt_fact = sqrt(factorial(a_val) * factorial(b_val));
    phase_arg = -b_val * rho;

    % Inner sum over l = 0 .. (n - j)
    inner = zeros(size(g));
    for l = 0:(n - j)
        coeff_l = ((-1)^l) / (factorial(a_val - l) * factorial(b_val - l) * factorial(l));
        inner = inner + coeff_l ./ eps_safe.^l;
    end

    term = sgn * c_i * c_f * eps_pow .* exp_eps * sqrt_fact .* inner;
    total_re = total_re + term .* cos(phase_arg);
    total_im = total_im + term .* sin(phase_arg);
end

P = total_re.^2 + total_im.^2;   % |total|^2
P(eps <= 0) = 0;
P = min(max(P, 0), 1);

end


% =====================================================================
% Zelechow orthogonal transformation matrix C^(n)
% =====================================================================
function C = zelechow_C_matrix(n_total)
% ZELECHOW_C_MATRIX  Orthogonal transformation matrix C^(n).
%
%   C(k, q) with 1-based indexing corresponds to
%   C^(n)_{k,q} in the Zelechow convention.

dim = n_total + 1;
C = zeros(dim, dim);
for kk = 1:dim
    for qq = 1:dim
        pref = 2^(-n_total / 2);
        pref = pref * sqrt(nchoosek_safe(n_total, kk-1) * nchoosek_safe(n_total, qq-1));
        s = 0;
        for nu = 0:(qq-1)
            s = s + ((-1)^(n_total + 1 - nu)) ...
                * nchoosek_safe(n_total - kk + 1, qq - nu - 1) ...
                * nchoosek_safe(kk - 1, nu);
        end
        C(kk, qq) = pref * s;
    end
end
end


% =====================================================================
% Safe binomial coefficient (returns 0 for invalid inputs)
% =====================================================================
function val = nchoosek_safe(n, k)
if k < 0 || k > n || n < 0
    val = 0;
else
    val = nchoosek(n, k);
end
end


% =====================================================================
% Numerically-safe cosh^2(a) / sinh^2(b)
% =====================================================================
function result = safe_cosh2_over_sinh2(a, b)

a_abs = abs(a);
b_abs = abs(b);
result = zeros(size(a));

both_large = (a_abs > 50) & (b_abs > 50);
result(both_large) = exp(2 * (a_abs(both_large) - b_abs(both_large)));

a_large = (a_abs > 50) & (b_abs <= 50);
if any(a_large)
    c2a = 0.25 * exp(2 * a_abs(a_large));
    s2b = sinh(b(a_large)).^2;
    s2b = max(s2b, 1e-100);
    result(a_large) = c2a ./ s2b;
end

b_large = (a_abs <= 50) & (b_abs > 50);
if any(b_large)
    c2a = cosh(a(b_large)).^2;
    s2b = 0.25 * exp(2 * b_abs(b_large));
    result(b_large) = c2a ./ s2b;
end

both_small = (a_abs <= 50) & (b_abs <= 50);
if any(both_small)
    s2b = sinh(b(both_small)).^2;
    s2b = max(s2b, 1e-100);
    result(both_small) = cosh(a(both_small)).^2 ./ s2b;
end

end

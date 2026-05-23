function [R_donor, R_partner_E] = ch4_vv_partner_source(T, x, AD, n_scale, n_steps, donor_delta, partner_mode, steric_vv, gamma_vv)
%CH4_VV_PARTNER_SOURCE Pair-resolved source for VV channels with a partner excitation.
%   donor_delta maps the donor state i -> f; partner_mode is 2 or 4 and
%   the partner transition is k_q -> k_q + 1.

if nargin < 5 || isempty(n_steps)
    n_steps = 6000;
end
if nargin < 9 || isempty(gamma_vv)
    gamma_vv = 0.5;
end

N = numel(AD.e1234);
R_donor = zeros(N, 1);
R_partner_E = 0;

if partner_mode ~= 2 && partner_mode ~= 4
    error('ch4_vv_partner_source:bad_partner_mode', ...
        'partner_mode must be 2 or 4.');
end

inds = AD.inds;
map = containers.Map('KeyType', 'char', 'ValueType', 'double');
for r = 1:N
    key = sprintf('%d,%d,%d,%d', inds(r, 1), inds(r, 2), inds(r, 3), inds(r, 4));
    map(key) = r;
end

qmax = AD.lch4(partner_mode) - 1;
donor_active = find(donor_delta ~= 0);   % VV rate depends only on active+partner modes
partner_low_pop = zeros(qmax + 1, 1);
partner_high_pop = zeros(qmax + 1, 1);
for p = 1:N
    q = inds(p, partner_mode);
    if q <= qmax - 1
        partner_low_pop(q + 1) = partner_low_pop(q + 1) + x(p);
    end
    if q >= 1
        partner_high_pop(q + 1) = partner_high_pop(q + 1) + x(p);
    end
end

eps = mode_energy(partner_mode, AD);
rate_cache = containers.Map('KeyType', 'char', 'ValueType', 'double');

for r = 1:N
    si = inds(r, :);
    sf = si + donor_delta;
    if any(sf < 0)
        continue
    end
    keyf = sprintf('%d,%d,%d,%d', sf(1), sf(2), sf(3), sf(4));
    if ~isKey(map, keyf)
        continue
    end
    dst = map(keyf);

    for qp = 0:(qmax - 1)
        x_low = partner_low_pop(qp + 1);
        x_high = partner_high_pop(qp + 2);
        if x_low == 0 && x_high == 0
            continue
        end

        sk = [0, 0, 0, 0];
        skf = [0, 0, 0, 0];
        sk(partner_mode) = qp;
        skf(partner_mode) = qp + 1;

        key = sprintf('%s%d_%d', sprintf('%d_', si(donor_active)), partner_mode, qp);
        if isKey(rate_cache, key)
            kf = rate_cache(key);
        else
            kf = rates_fho_vv(T, si, sf, sk, skf, steric_vv, ...
                AD.fho_alpha, AD.fho_e_m, 'n_steps', n_steps, 'gamma', gamma_vv);
            rate_cache(key) = kf;
        end

        kr = ch4_vv_pair_backward(kf, r, dst, sk, skf, T, AD);
        flux = n_scale * (x(r) * x_low * kf - x(dst) * x_high * kr);

        R_donor(r) = R_donor(r) - flux;
        R_donor(dst) = R_donor(dst) + flux;
        R_partner_E = R_partner_E + eps * flux;
    end
end

end

function kb = ch4_vv_pair_backward(kf, i, f, k, kf_state, T, AD)
s_i = AD.stw(i);
s_f = AD.stw(f);
s_k = ch4_state_weight(k);
s_kf = ch4_state_weight(kf_state);

delta_e = AD.e1234(f) - AD.e1234(i) + ch4_state_energy_delta(k, kf_state, AD);
kb = kf * (s_i * s_k) / (s_f * s_kf) * exp(delta_e / (AD.k * T));
end

function s = ch4_state_weight(state)
s = (state(2) + 1) * (state(3) + 1) * (state(3) + 2) * ...
    (state(4) + 1) * (state(4) + 2) / 4;
end

function de = ch4_state_energy_delta(state_i, state_f, AD)
de = 0;
for mode = 1:4
    de = de + (state_f(mode) - state_i(mode)) * mode_energy(mode, AD);
end
end

function eps = mode_energy(mode, AD)
switch mode
    case 1
        eps = AD.e1000 - AD.e0000;
    case 2
        eps = AD.e0100 - AD.e0000;
    case 3
        eps = AD.e0010 - AD.e0000;
    case 4
        eps = AD.e0001 - AD.e0000;
    otherwise
        error('ch4_vv_partner_source:bad_mode', 'mode must be 1, 2, 3, or 4.');
end
end

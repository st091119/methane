function [R_donor, R_partner_E] = ch4_vv34_4d_source(T, x, AD, n_scale, n_steps)
%CH4_VV34_4D_SOURCE Source for CH4(i)+CH4(k4)->CH4(i3-1,i4+1)+CH4(k4+1).
%   R_donor is the state source for the donor molecule only. R_partner_E is
%   the matching vibrational-energy source of the partner molecule.

if nargin < 5
    n_steps = 6000;
end

if ~isfield(AD, 'indvv34') || ~isfield(AD, 'indl_vt4')
    AD = ch4_relax_topology(AD);
end

N = numel(AD.e1234);
R_donor = zeros(N, 1);
R_partner_E = 0;

eps4 = AD.e0001 - AD.e0000;
Lmx = AD.lch4(4) - 1;

partner_low_pop = zeros(Lmx + 1, 1);
partner_high_pop = zeros(Lmx + 1, 1);
partner_src = find(~isnan(AD.indl_vt4));
for ii = 1:numel(partner_src)
    p = partner_src(ii);
    pf = AD.indl_vt4(p);
    Lp = AD.inds(p, 4);
    Lpf = AD.inds(pf, 4);
    partner_low_pop(Lp + 1) = partner_low_pop(Lp + 1) + x(p);
    partner_high_pop(Lpf + 1) = partner_high_pop(Lpf + 1) + x(pf);
end

src = find(~isnan(AD.indvv34));
rate_cache = containers.Map('KeyType', 'char', 'ValueType', 'double');

for ii = 1:numel(src)
    r = src(ii);
    dst = AD.indvv34(r);
    si = AD.inds(r, :);
    sf = AD.inds(dst, :);

    for Lp = 0:(Lmx - 1)
        x_low = partner_low_pop(Lp + 1);
        x_high = partner_high_pop(Lp + 2);
        if x_low == 0 && x_high == 0
            continue
        end

        sk = [0, 0, 0, Lp];
        skf = [0, 0, 0, Lp + 1];
        key = sprintf('%d_%d_%d', si(3), si(4), Lp);
        if isKey(rate_cache, key)
            kf = rate_cache(key);
        else
            kf = rates_fho_vv(T, si, sf, sk, skf, ...
                AD.fho_steric_vv_34_4d, ...
                AD.fho_alpha, AD.fho_e_m, 'n_steps', n_steps);
            rate_cache(key) = kf;
        end

        kr = ch4_vv_pair_backward(kf, r, dst, sk, skf, T, AD);
        flux = n_scale * (x(r) * x_low * kf - x(dst) * x_high * kr);

        R_donor(r) = R_donor(r) - flux;
        R_donor(dst) = R_donor(dst) + flux;
        R_partner_E = R_partner_E + eps4 * flux;
    end
end

end

function kb = ch4_vv_pair_backward(kf, i, f, k, kf_state, T, AD)
% Backward pair rate for i+k -> f+kf_state.
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
eps1 = AD.e1000 - AD.e0000;
eps2 = AD.e0100 - AD.e0000;
eps3 = AD.e0010 - AD.e0000;
eps4 = AD.e0001 - AD.e0000;
dq = state_f - state_i;
de = dq(1) * eps1 + dq(2) * eps2 + dq(3) * eps3 + dq(4) * eps4;
end

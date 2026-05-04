function B = ch4_B_vv34_4d(T, T3star, T4star, AD, n_steps)
%CH4_B_VV34_4D  Macroscopic B_{34,4}^d(T, T3*, T4*) [m^3/s].
%   Sum / (Z3*Z4^2), forward rate minus reverse, same FHO states as old pair model.

if nargin < 5
    n_steps = 6000;
end

if ~exist('rates_fho_vv', 'file')
    addpath(fullfile(fileparts(mfilename('fullpath')), 'fho_model'));
end

kB = AD.k;
eps3 = AD.e0010 - AD.e0000;
eps4 = AD.e0001 - AD.e0000;

T3s = max(T3star, 1e-12); 
T4s = max(T4star, 1e-12);

i3_max = AD.lch4(3) - 1;
i4_max = AD.lch4(4) - 1;

Z3 = 0;
for i3 = 0:i3_max
    Z3 = Z3 + stat_w_methane(i3, 3) * exp(-i3 * eps3 / (kB * T3s));
end
Z4 = 0;
for i4 = 0:i4_max
    Z4 = Z4 + stat_w_methane(i4, 4) * exp(-i4 * eps4 / (kB * T4s));
end

den = Z3 * Z4^2;

cache_f = containers.Map('KeyType', 'char', 'ValueType', 'double');

steric = AD.fho_steric_vv_34_4d;
alpha = AD.fho_alpha;
e_m = AD.fho_e_m;

acc = 0;
for i3 = 0:i3_max
    for i4 = 0:i4_max
        for j4 = 0:i4_max
            wpop = stat_w_methane(i3, 3) * stat_w_methane(i4, 4) * stat_w_methane(j4, 4);
            boltz = exp(-i3 * eps3 / (kB * T3s) - (i4 + j4) * eps4 / (kB * T4s));
            w = wpop * boltz;

            kf = 0;
            if i3 >= 1 && i4 < i4_max && j4 < i4_max 
                keyf = sprintf('%d_%d_%d', i3, i4, j4);
                if isKey(cache_f, keyf)
                    kf = cache_f(keyf);
                else
                    si = [0, 0, i3, i4];
                    sf = [0, 0, i3 - 1, i4 + 1];
                    sk = [0, 0, 0, j4];
                    skf = [0, 0, 0, j4 + 1];
                    kf = rates_fho_vv(T, si, sf, sk, skf, steric, alpha, e_m, 'n_steps', n_steps);
                    cache_f(keyf) = kf;
                end
            end

            kr = 0;
            if i3 < i3_max && i4 >= 1 && j4 >= 1
                si_rev = [0, 0, i3 + 1, i4 - 1];
                sf_rev = [0, 0, i3, i4];
                sk_rev = [0, 0, 0, j4 - 1];
                skf_rev = [0, 0, 0, j4];
                keyr = sprintf('%d_%d_%d', si_rev(3), si_rev(4), sk_rev(4));
                if isKey(cache_f, keyr)
                    kf_rev = cache_f(keyr);
                else
                    kf_rev = rates_fho_vv(T, si_rev, sf_rev, sk_rev, skf_rev, steric, alpha, e_m, 'n_steps', n_steps);
                    cache_f(keyr) = kf_rev;
                end
                kr = ch4_pair_backward(kf_rev, si_rev, sf_rev, sk_rev, skf_rev, T, AD);
            end

            acc = acc + w * (kf - kr);
        end
    end
end

B = acc / den;
end

function sw = stat_w_methane(level, mode)
switch mode
    case 3
        sw = 0.5 * (level + 1) * (level + 2);
    case 4
        sw = 0.5 * (level + 1) * (level + 2);
    otherwise
        error('ch4_B_vv34_4d:invalid_mode', 'mode must be 3 or 4.');
end
end

function kb = ch4_pair_backward(kf, state_i, state_f, state_k, state_kf, T, AD)
s_i = stat_w_pair(state_i);
s_f = stat_w_pair(state_f);
s_k = stat_w_pair(state_k);
s_kf = stat_w_pair(state_kf);
de = ch4_state_energy(state_f, AD) - ch4_state_energy(state_i, AD) + ...
    ch4_state_energy(state_kf, AD) - ch4_state_energy(state_k, AD);
kb = kf * (s_i * s_k) / (s_f * s_kf) * exp(de / (AD.k * T));
end

function s = stat_w_pair(state)
s = (state(2) + 1) * (state(3) + 1) * (state(3) + 2) * ...
    (state(4) + 1) * (state(4) + 2) / 4;
end

function e = ch4_state_energy(state, AD)
eps1 = AD.e1000 - AD.e0000;
eps2 = AD.e0100 - AD.e0000;
eps3 = AD.e0010 - AD.e0000;
eps4 = AD.e0001 - AD.e0000;
e = state(1) * eps1 + state(2) * eps2 + state(3) * eps3 + state(4) * eps4;
end

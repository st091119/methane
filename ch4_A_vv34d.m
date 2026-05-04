function A = ch4_A_vv34d(T, T3star, T4star, AD, n_steps)
%CH4_A_VV34D  Macroscopic A_{34}^d(T, T3*, T4*) [m^3/s].
%   A = (1/(Z3*Z4)) * sum_{i3,i4} s3*s4*exp(-i3*eps3/(k*T3*)-i4*eps4/(k*T4*))*(k_f-k_r)
%   Forward: (i3,i4)->(i3-1,i4+2), reverse: (i3,i4)->(i3+1,i4-2).

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
    Z3 = Z3 + stat_w_methane(i3) * exp(-i3 * eps3 / (kB * T3s));
end
Z4 = 0;
for i4 = 0:i4_max
    Z4 = Z4 + stat_w_methane(i4) * exp(-i4 * eps4 / (kB * T4s));
end

den = Z3 * Z4;
if den <= 0 || ~isfinite(den)
    A = 0;
    return
end

cache_f = containers.Map('KeyType', 'char', 'ValueType', 'double');

steric = AD.fho_steric_vv_34d;
alpha = AD.fho_alpha;
e_m = AD.fho_e_m;
gamma = AD.fho_gamma_vv34d;

acc = 0;
sk0 = [0, 0, 0, 0];
for i3 = 0:i3_max
    for i4 = 0:i4_max
        w = stat_w_methane(i3) * stat_w_methane(i4) * ...
            exp(-i3 * eps3 / (kB * T3s) - i4 * eps4 / (kB * T4s));

        kf = 0;
        if i3 >= 1 && i4 <= i4_max - 2
            keyf = sprintf('%d_%d', i3, i4);
            if isKey(cache_f, keyf)
                kf = cache_f(keyf);
            else
                si = [0, 0, i3, i4];
                sf = [0, 0, i3 - 1, i4 + 2];
                kf = rates_fho_vv(T, si, sf, sk0, sk0, steric, alpha, e_m, ...
                    'n_steps', n_steps, 'gamma', gamma);
                cache_f(keyf) = kf;
            end
        end

        kr = 0;
        if i3 <= i3_max - 1 && i4 >= 2
            si_rev = [0, 0, i3 + 1, i4 - 2];
            sf_rev = [0, 0, i3, i4];
            keyr = sprintf('%d_%d', si_rev(3), si_rev(4));
            if isKey(cache_f, keyr)
                kf_rev = cache_f(keyr);
            else
                kf_rev = rates_fho_vv(T, si_rev, sf_rev, sk0, sk0, steric, alpha, e_m, ...
                    'n_steps', n_steps, 'gamma', gamma);
                cache_f(keyr) = kf_rev;
            end
            kr = ch4_state_backward(kf_rev, si_rev, sf_rev, T, AD);
        end

        acc = acc + w * (kf - kr);
    end
end

A = acc / den;
end

function sw = stat_w_methane(level)
sw = 0.5 * (level + 1) * (level + 2);
end

function kb = ch4_state_backward(kf, state_i, state_f, T, AD)
s_i = stat_w_methane(state_i(3)) * stat_w_methane(state_i(4));
s_f = stat_w_methane(state_f(3)) * stat_w_methane(state_f(4));
de = ch4_state_energy(state_f, AD) - ch4_state_energy(state_i, AD);
kb = kf * (s_i / s_f) * exp(de / (AD.k * T));
end

function e = ch4_state_energy(state, AD)
eps3 = AD.e0010 - AD.e0000;
eps4 = AD.e0001 - AD.e0000;
e = state(3) * eps3 + state(4) * eps4;
end

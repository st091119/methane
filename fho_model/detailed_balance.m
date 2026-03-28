function varargout = detailed_balance(func_name, varargin)
% DETAILED_BALANCE  Detailed balance calculations for VT and VV processes.
%
%   Collection of functions for calculating backward rate coefficients
%   from forward rates using the detailed balance principle.
%
%   Available functions (call via detailed_balance('func_name', args...)):
%
%     s = detailed_balance('stat_weight', state)
%         Statistical weight s for CH4 state [i1, i2, i3, i4].
%         s = (i2+1)(i3+1)(i3+2)(i4+1)(i4+2) / 4
%
%     eps = detailed_balance('energy', state, AD)
%         Vibrational energy [J] of CH4 state from ground state (harmonic).
%
%     db = detailed_balance('db_factor_vt', state_i, state_f, T, AD)
%         Detailed balance ratio k_{f->i} / k_{i->f} for VT process.
%
%     db = detailed_balance('db_factor_vv', state_i, state_f, state_k, state_kf, T, AD)
%         Detailed balance ratio k_{f,kf->i,k} / k_{i,k->f,kf} for VV process.
%
%     k_bwd = detailed_balance('rate_vt_backward', T, state_i, state_f, k_fwd, AD)
%         Backward VT rate from forward rate.
%
%     k_bwd = detailed_balance('rate_vv_backward', T, si, sf, sk, skf, k_fwd, AD)
%         Backward VV rate from forward rate.
%
%   Requires: input_data.m

switch func_name
    case 'stat_weight'
        varargout{1} = ch4_stat_weight(varargin{:});
    case 'energy'
        varargout{1} = ch4_energy(varargin{:});
    case 'db_factor_vt'
        varargout{1} = db_factor_vt(varargin{:});
    case 'db_factor_vv'
        varargout{1} = db_factor_vv(varargin{:});
    case 'rate_vt_backward'
        varargout{1} = rate_vt_backward(varargin{:});
    case 'rate_vv_backward'
        varargout{1} = rate_vv_backward(varargin{:});
    otherwise
        error('Unknown function: %s', func_name);
end

end


% =====================================================================
% Statistical weight for CH4 state
% =====================================================================
function s = ch4_stat_weight(state)
% CH4_STAT_WEIGHT  Statistical weight s for CH4 state [i1, i2, i3, i4].
%   s = (i2+1)(i3+1)(i3+2)(i4+1)(i4+2) / 4

i2 = state(2);
i3 = state(3);
i4 = state(4);
s = (i2 + 1) * (i3 + 1) * (i3 + 2) * (i4 + 1) * (i4 + 2) / 4.0;

end


% =====================================================================
% Vibrational energy from ground state
% =====================================================================
function eps = ch4_energy(state, AD)
% CH4_ENERGY  Vibrational energy [J] of CH4 state from ground state.
%   eps = hc * sum_m (omega_m * i_m)  [harmonic oscillator]

hc    = AD.h * AD.c;   % [J*cm]
omega = AD.omega;       % [cm^-1]

eps = hc * (omega(1) * state(1) + ...
            omega(2) * state(2) + ...
            omega(3) * state(3) + ...
            omega(4) * state(4));

end


% =====================================================================
% VT detailed balance factor
% =====================================================================
function db = db_factor_vt(state_i, state_f, T, AD)
% DB_FACTOR_VT  Detailed balance ratio k_{f->i} / k_{i->f} for VT.
%
%   Forward:  CH4(state_i) + M -> CH4(state_f) + M
%   Returns:  k_{f->i}(T) / k_{i->f}(T)
%
%   db = (s_i / s_f) * exp((eps_f - eps_i) / (k_B * T))

s_i   = ch4_stat_weight(state_i);
s_f   = ch4_stat_weight(state_f);
eps_i = ch4_energy(state_i, AD);
eps_f = ch4_energy(state_f, AD);

db = (s_i / s_f) * exp((eps_f - eps_i) / (AD.k * T));

end


% =====================================================================
% VV detailed balance factor
% =====================================================================
function db = db_factor_vv(state_i, state_f, state_k, state_kf, T, AD)
% DB_FACTOR_VV  Detailed balance ratio for VV process.
%
%   Forward:  CH4(i) + CH4(k) -> CH4(f) + CH4(kf)
%   Returns:  k_{f,kf->i,k}(T) / k_{i,k->f,kf}(T)
%
%   db = (s_i * s_k) / (s_f * s_kf) * exp((eps_f + eps_kf - eps_i - eps_k) / (k_B * T))

s_i   = ch4_stat_weight(state_i);
s_f   = ch4_stat_weight(state_f);
s_k   = ch4_stat_weight(state_k);
s_kf  = ch4_stat_weight(state_kf);

eps_i  = ch4_energy(state_i,  AD);
eps_f  = ch4_energy(state_f,  AD);
eps_k  = ch4_energy(state_k,  AD);
eps_kf = ch4_energy(state_kf, AD);

db = (s_i * s_k) / (s_f * s_kf) * ...
     exp((eps_f + eps_kf - eps_i - eps_k) / (AD.k * T));

end


% =====================================================================
% Backward VT rate from forward rate
% =====================================================================
function k_bwd = rate_vt_backward(T, state_i, state_f, k_fwd, AD)
% RATE_VT_BACKWARD  Backward VT rate via detailed balance.
%
%   Forward:  CH4(state_i) + M -> CH4(state_f) + M  with rate k_fwd
%   Returns:  k_{f->i}(T)  [same units as k_fwd]

k_bwd = k_fwd * db_factor_vt(state_i, state_f, T, AD);

end


% =====================================================================
% Backward VV rate from forward rate
% =====================================================================
function k_bwd = rate_vv_backward(T, state_i, state_f, state_k, state_kf, k_fwd, AD)
% RATE_VV_BACKWARD  Backward VV rate via detailed balance.
%
%   Forward:  CH4(i) + CH4(k) -> CH4(f) + CH4(kf)  with rate k_fwd
%   Returns:  k_{f,kf->i,k}(T)  [same units as k_fwd]

k_bwd = k_fwd * db_factor_vv(state_i, state_f, state_k, state_kf, T, AD);

end

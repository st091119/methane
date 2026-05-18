function terms = ch4_relax_terms_simplified(T, varargin)
%CH4_RELAX_TERMS_SIMPLIFIED Reduced CH4 relaxation terms for 2T and 3T models.
%
%   terms = ch4_relax_terms_simplified(T, Tv, AD)
%   terms = ch4_relax_terms_simplified(T, T13, T24, AD)
%
% Returned relaxation terms are nondimensional RHS moments, matching the
% scaling used in rpart_mt_sts.m and rpart_mt_3t_sts.m:
%   R_nd = R_dim * tau / (n0*k*T0).
%
% VT2 and VT4 use the regression models from regression.tex:
%   S_m = (1/Tstar - 1/T) * exp(P_m(1000/T, 1000/Tstar)).
%
% VV34s, VV34d, and VV34_4d use the reduced sums from formulas_scheme1.tex.
%
% Cache helpers:
%   ch4_relax_terms_simplified('clear_cache')
%   coeffs = ch4_relax_terms_simplified('coefficients')
%
% Options:
%   'processes'  cell/string list, default {'VT2','VT4','VV34s','VV34d','VV34_4d'}
%   'n_steps_vv' number of FHO quadrature steps for VV rates
%   'use_cache'  true/false

if ischar(T) || isstring(T)
    cmd = lower(string(T));
    switch cmd
        case "clear_cache"
            cached_vv_rate('clear');
            terms = [];
            return
        case "coefficients"
            terms = regression_coefficients();
            return
        otherwise
            error('ch4_relax_terms_simplified:unknown_command', ...
                'Unknown command: %s', cmd);
    end
end

if numel(varargin) < 2
    error('ch4_relax_terms_simplified:not_enough_inputs', ...
        'Use (T,Tv,AD) for 2T or (T,T13,T24,AD) for 3T.');
end

if isstruct(varargin{2})
    model = '2T';
    Tv = varargin{1};
    AD = varargin{2};
    opt_args = varargin(3:end);
else
    if numel(varargin) < 3 || ~isstruct(varargin{3})
        error('ch4_relax_terms_simplified:bad_inputs', ...
            'Use (T,Tv,AD) for 2T or (T,T13,T24,AD) for 3T.');
    end
    model = '3T';
    T13 = varargin{1};
    T24 = varargin{2};
    AD = varargin{3};
    opt_args = varargin(4:end);
end

opts = parse_options(opt_args{:});
AD = ensure_relax_defaults(AD);

eps = mode_energies(AD);
scale_nd = AD.n0 * AD.tau / (AD.k * AD.T0);
scale_dim = AD.n0^2;

terms = struct();
terms.model = model;
terms.T = T;
terms.n0 = AD.n0;
terms.tau = AD.tau;
terms.T0 = AD.T0;
terms.scale_nd = scale_nd;
terms.scale_dim = scale_dim;
terms.units = 'nondimensional RHS moment';

switch model
    case '2T'
        terms.Tv = Tv;
        terms = fill_2t_terms(terms, T, Tv, AD, eps, scale_nd, scale_dim, opts);
    case '3T'
        terms.T13 = T13;
        terms.T24 = T24;
        terms = fill_3t_terms(terms, T, T13, T24, AD, eps, scale_nd, scale_dim, opts);
end
end

function opts = parse_options(varargin)
opts = struct();
opts.n_steps_vv = 6000;
opts.use_cache = true;
opts.processes = {'vt2', 'vt4', 'vv34s', 'vv34d', 'vv34_4d'};

if mod(numel(varargin), 2) ~= 0
    error('ch4_relax_terms_simplified:bad_options', ...
        'Options must be name-value pairs.');
end

for i = 1:2:numel(varargin)
    name = lower(string(varargin{i}));
    value = varargin{i + 1};
    switch name
        case "n_steps_vv"
            opts.n_steps_vv = value;
        case "use_cache"
            opts.use_cache = logical(value);
        case "processes"
            opts.processes = canonical_processes(value);
        otherwise
            error('ch4_relax_terms_simplified:unknown_option', ...
                'Unknown option: %s', name);
    end
end
end

function processes = canonical_processes(value)
if ischar(value) || isstring(value)
    value = cellstr(string(value));
elseif ~iscell(value)
    error('ch4_relax_terms_simplified:bad_processes', ...
        'Processes must be a string or a cell array of strings.');
end
if isscalar(value) && iscell(value{1})
    value = value{1};
end

processes = cell(size(value));
for i = 1:numel(value)
    name = lower(char(string(value{i})));
    name = strrep(name, '^', '');
    name = strrep(name, '{', '');
    name = strrep(name, '}', '');
    name = strrep(name, '-', '');
    name = strrep(name, ',', '');
    name = strrep(name, ' ', '');
    switch name
        case {'vt2'}
            processes{i} = 'vt2';
        case {'vt4'}
            processes{i} = 'vt4';
        case {'vv34s', 'vvs34', 'vvs_34', 'vv_s_34'}
            processes{i} = 'vv34s';
        case {'vv34d', 'vvd34', 'vvd_34'}
            processes{i} = 'vv34d';
        case {'vv34_4d', 'vv344d', 'vvd344', 'vvd344d', 'vvd34_4', 'vvd34_4d'}
            processes{i} = 'vv34_4d';
        otherwise
            error('ch4_relax_terms_simplified:unknown_process', ...
                'Unknown relaxation process: %s', char(string(value{i})));
    end
end
processes = unique(processes, 'stable');
end

function tf = has_process(opts, name)
tf = any(strcmpi(opts.processes, name));
end

function AD = ensure_relax_defaults(AD)
if ~exist('rates_fho_vv', 'file')
    addpath(fullfile(fileparts(mfilename('fullpath')), 'fho_model'));
end

if ~isfield(AD, 'T0')
    AD.T0 = 1000;
end
if ~isfield(AD, 'p0')
    AD.p0 = 101325;
end
if ~isfield(AD, 'n0')
    AD.n0 = AD.p0 / (AD.k * AD.T0);
end
if ~isfield(AD, 'tau')
    sigma0 = pi * AD.r0^2;
    AD.tau = (4 * AD.n0 * sigma0 * sqrt(AD.k * AD.T0 / (pi * AD.m)))^(-1);
end

if ~isfield(AD, 'fho_alpha')
    AD.fho_alpha = 5.174e10;
end
if ~isfield(AD, 'fho_e_m')
    AD.fho_e_m = 776.4;
end
if ~isfield(AD, 'fho_steric_vv_34s')
    AD.fho_steric_vv_34s = 1.0;
end
if ~isfield(AD, 'fho_steric_vv_34d')
    AD.fho_steric_vv_34d = 0.06;
end
if ~isfield(AD, 'fho_gamma_vv34d')
    AD.fho_gamma_vv34d = 1.0;
end
if ~isfield(AD, 'fho_steric_vv_34_4d')
    AD.fho_steric_vv_34_4d = 1.0;
end
if ~isfield(AD, 'fho_steric_vv_32_2d')
    AD.fho_steric_vv_32_2d = 1.0;
end
if ~isfield(AD, 'fho_steric_vv_34_2d')
    AD.fho_steric_vv_34_2d = 1.0;
end
if ~isfield(AD, 'fho_steric_vv_32_4d')
    AD.fho_steric_vv_32_4d = 1.0;
end
end

function eps = mode_energies(AD)
eps.eps1 = AD.e1000 - AD.e0000;
eps.eps2 = AD.e0100 - AD.e0000;
eps.eps3 = AD.e0010 - AD.e0000;
eps.eps4 = AD.e0001 - AD.e0000;
end

function terms = fill_2t_terms(terms, T, Tv, AD, eps, scale_nd, scale_dim, opts)
terms.VT2 = source_entry(0, 0, 0);
terms.VT4 = source_entry(0, 0, 0);
terms.VV34s = source_entry(0, 0, 0);
terms.VV34d = source_entry(0, 0, 0);
terms.VV34_4d = source_entry(0, 0, 0);

if has_process(opts, 'vt2')
    S2 = vt_regression(2, T, Tv);
    terms.VT2 = source_entry(scale_nd * eps.eps2 * S2, scale_dim * eps.eps2 * S2, S2);
end
if has_process(opts, 'vt4')
    S4 = vt_regression(4, T, Tv);
    terms.VT4 = source_entry(scale_nd * eps.eps4 * S4, scale_dim * eps.eps4 * S4, S4);
end
if has_process(opts, 'vv34s')
    A34s = vv34s_sum(T, Tv, Tv, AD, opts);
    terms.VV34s = source_entry(scale_nd * (eps.eps4 - eps.eps3) * A34s, ...
        scale_dim * (eps.eps4 - eps.eps3) * A34s, A34s);
end
if has_process(opts, 'vv34d')
    A34d = vv34d_sum(T, Tv, Tv, AD, opts);
    terms.VV34d = source_entry(scale_nd * (2 * eps.eps4 - eps.eps3) * A34d, ...
        scale_dim * (2 * eps.eps4 - eps.eps3) * A34d, A34d);
end
if has_process(opts, 'vv34_4d')
    B34d4 = vv34_4d_sum(T, Tv, Tv, AD, opts);
    terms.VV34_4d = source_entry(scale_nd * (2 * eps.eps4 - eps.eps3) * B34d4, ...
        scale_dim * (2 * eps.eps4 - eps.eps3) * B34d4, B34d4);
end

terms.Rvibr = terms.VT2.nd + terms.VT4.nd + terms.VV34s.nd + ...
    terms.VV34d.nd + terms.VV34_4d.nd;
terms.Rvibr_dim = terms.VT2.dim + terms.VT4.dim + terms.VV34s.dim + ...
    terms.VV34d.dim + terms.VV34_4d.dim;
end

function terms = fill_3t_terms(terms, T, T13, T24, AD, eps, scale_nd, scale_dim, opts)
terms.VT2.R13 = 0;
terms.VT2.R24 = 0;
terms.VT2.R13_dim = 0;
terms.VT2.R24_dim = 0;
terms.VT2.reduced_sum = 0;

terms.VT4.R13 = 0;
terms.VT4.R24 = 0;
terms.VT4.R13_dim = 0;
terms.VT4.R24_dim = 0;
terms.VT4.reduced_sum = 0;

terms.VV34s.R13 = 0;
terms.VV34s.R24 = 0;
terms.VV34s.R13_dim = 0;
terms.VV34s.R24_dim = 0;
terms.VV34s.reduced_sum = 0;

terms.VV34d.R13 = 0;
terms.VV34d.R24 = 0;
terms.VV34d.R13_dim = 0;
terms.VV34d.R24_dim = 0;
terms.VV34d.reduced_sum = 0;

terms.VV34_4d.R13 = 0;
terms.VV34_4d.R24 = 0;
terms.VV34_4d.R13_dim = 0;
terms.VV34_4d.R24_dim = 0;
terms.VV34_4d.reduced_sum = 0;

if has_process(opts, 'vt2')
    S2 = vt_regression(2, T, T24);
    terms.VT2.R24 = scale_nd * eps.eps2 * S2;
    terms.VT2.R24_dim = scale_dim * eps.eps2 * S2;
    terms.VT2.reduced_sum = S2;
end
if has_process(opts, 'vt4')
    S4 = vt_regression(4, T, T24);
    terms.VT4.R24 = scale_nd * eps.eps4 * S4;
    terms.VT4.R24_dim = scale_dim * eps.eps4 * S4;
    terms.VT4.reduced_sum = S4;
end
if has_process(opts, 'vv34s')
    A34s = vv34s_sum(T, T13, T24, AD, opts);
    terms.VV34s.R13 = -scale_nd * eps.eps3 * A34s;
    terms.VV34s.R24 = scale_nd * eps.eps4 * A34s;
    terms.VV34s.R13_dim = -scale_dim * eps.eps3 * A34s;
    terms.VV34s.R24_dim = scale_dim * eps.eps4 * A34s;
    terms.VV34s.reduced_sum = A34s;
end
if has_process(opts, 'vv34d')
    A34d = vv34d_sum(T, T13, T24, AD, opts);
    terms.VV34d.R13 = -scale_nd * eps.eps3 * A34d;
    terms.VV34d.R24 = 2 * scale_nd * eps.eps4 * A34d;
    terms.VV34d.R13_dim = -scale_dim * eps.eps3 * A34d;
    terms.VV34d.R24_dim = 2 * scale_dim * eps.eps4 * A34d;
    terms.VV34d.reduced_sum = A34d;
end
if has_process(opts, 'vv34_4d')
    B34d4 = vv34_4d_sum(T, T13, T24, AD, opts);
    terms.VV34_4d.R13 = -scale_nd * eps.eps3 * B34d4;
    terms.VV34_4d.R24 = 2 * scale_nd * eps.eps4 * B34d4;
    terms.VV34_4d.R13_dim = -scale_dim * eps.eps3 * B34d4;
    terms.VV34_4d.R24_dim = 2 * scale_dim * eps.eps4 * B34d4;
    terms.VV34_4d.reduced_sum = B34d4;
end

terms.R13 = terms.VT2.R13 + terms.VT4.R13 + terms.VV34s.R13 + ...
    terms.VV34d.R13 + terms.VV34_4d.R13;
terms.R24 = terms.VT2.R24 + terms.VT4.R24 + terms.VV34s.R24 + ...
    terms.VV34d.R24 + terms.VV34_4d.R24;
terms.R13_dim = terms.VT2.R13_dim + terms.VT4.R13_dim + terms.VV34s.R13_dim + ...
    terms.VV34d.R13_dim + terms.VV34_4d.R13_dim;
terms.R24_dim = terms.VT2.R24_dim + terms.VT4.R24_dim + terms.VV34s.R24_dim + ...
    terms.VV34d.R24_dim + terms.VV34_4d.R24_dim;
end

function entry = source_entry(nd, dim, reduced_sum)
entry = struct();
entry.nd = nd;
entry.dim = dim;
entry.reduced_sum = reduced_sum;
end

function S = vt_regression(mode, T, Tstar)
coeffs = regression_coefficients();
switch mode
    case 2
        c = coeffs.vt2;
    case 4
        c = coeffs.vt4;
    otherwise
        error('ch4_relax_terms_simplified:bad_vt_mode', 'mode must be 2 or 4.');
end

theta = 1000;
u = theta / max(T, realmin);
v = theta / max(Tstar, realmin);
P = 0;
k = 1;
for d = 0:6
    for p = 0:d
        P = P + c(k) * u^p * v^(d - p);
        k = k + 1;
    end
end

drive = 1 / Tstar - 1 / T;
S = drive * exp(P);
end

function coeffs = regression_coefficients()
coeffs = struct();
coeffs.vt2 = [
    -3.304841133891e1
    -7.217904191638e0
    -1.043804734520e1
     5.343641347767e0
    -4.331466360961e-1
     5.805167325293e0
    -2.310393374327e0
    -3.311821797291e-3
    -2.881322842702e-3
    -2.219173274720e0
     5.677914344879e-1
     3.985332999075e-2
    -5.384829979249e-2
     3.941549367346e-2
     5.195352411197e-1
    -7.406972721173e-2
    -4.506694580732e-3
     2.511215095721e-3
     2.008857725780e-3
    -4.165168505388e-3
    -6.685027405438e-2
     3.996933748241e-3
    -3.285376853781e-4
     1.588899572278e-3
    -2.392888124160e-3
     1.632065925547e-3
    -3.681829015292e-4
     3.620567999665e-3
];

coeffs.vt4 = [
    -3.095247519915e1
    -9.630023306734e0
    -9.059542431293e0
     7.235892236171e0
    -4.584786974762e-1
     5.000509396068e0
    -3.131523109067e0
     5.011840743599e-2
     5.415716058063e-2
    -1.920389337523e0
     7.700721788399e-1
     1.187803498117e-2
    -4.446005645181e-2
     1.014208857510e-2
     4.529985123236e-1
    -1.001714694352e-1
    -9.724891760747e-4
     2.784366851462e-3
     2.606138045549e-3
    -5.560898731667e-4
    -5.852376786973e-2
     5.350160668027e-3
    -2.834986644878e-4
     8.655740294679e-4
    -1.469836649865e-3
     8.810098629606e-4
    -3.182132011700e-4
     3.160100934634e-3
];
end

function A = vv34s_sum(T, T3star, T4star, AD, opts)
eps = mode_energies(AD);
[Z3, Z4] = partition_34(T3star, T4star, eps, AD);
den = Z3 * Z4;
if den <= 0 || ~isfinite(den)
    A = 0;
    return
end

i3_max = AD.lch4(3) - 1;
i4_max = AD.lch4(4) - 1;
acc = 0;
for i3 = 0:i3_max
    for i4 = 0:i4_max
        w = stat_w3(i3) * stat_w4(i4) * ...
            exp(-i3 * eps.eps3 / (AD.k * T3star) - i4 * eps.eps4 / (AD.k * T4star));
        kp = 0;
        if i3 >= 1 && i4 <= i4_max - 1
            kp = cached_vv_rate('34s', T, i3, i4, 0, AD, opts);
        end
        km = 0;
        if i3 <= i3_max - 1 && i4 >= 1
            kprev = cached_vv_rate('34s', T, i3 + 1, i4 - 1, 0, AD, opts);
            km = detailed_balance_vv(kprev, i3 + 1, i4 - 1, 0, i3, i4, 0, T, eps, AD, false);
        end
        acc = acc + w * (kp - km);
    end
end
A = acc / den;
end

function A = vv34d_sum(T, T3star, T4star, AD, opts)
eps = mode_energies(AD);
[Z3, Z4] = partition_34(T3star, T4star, eps, AD);
den = Z3 * Z4;
if den <= 0 || ~isfinite(den)
    A = 0;
    return
end

i3_max = AD.lch4(3) - 1;
i4_max = AD.lch4(4) - 1;
acc = 0;
for i3 = 0:i3_max
    for i4 = 0:i4_max
        w = stat_w3(i3) * stat_w4(i4) * ...
            exp(-i3 * eps.eps3 / (AD.k * T3star) - i4 * eps.eps4 / (AD.k * T4star));
        kp = 0;
        if i3 >= 1 && i4 <= i4_max - 2
            kp = cached_vv_rate('34d', T, i3, i4, 0, AD, opts);
        end
        km = 0;
        if i3 <= i3_max - 1 && i4 >= 2
            kprev = cached_vv_rate('34d', T, i3 + 1, i4 - 2, 0, AD, opts);
            km = detailed_balance_vv(kprev, i3 + 1, i4 - 2, 0, i3, i4, 0, T, eps, AD, false);
        end
        acc = acc + w * (kp - km);
    end
end
A = acc / den;
end

function B = vv34_4d_sum(T, T3star, T4star, AD, opts)
eps = mode_energies(AD);
[Z3, Z4] = partition_34(T3star, T4star, eps, AD);
den = Z3 * Z4^2;
if den <= 0 || ~isfinite(den)
    B = 0;
    return
end

i3_max = AD.lch4(3) - 1;
i4_max = AD.lch4(4) - 1;
acc = 0;
for i3 = 0:i3_max
    for i4 = 0:i4_max
        for j4 = 0:i4_max
            w = stat_w3(i3) * stat_w4(i4) * stat_w4(j4) * ...
                exp(-i3 * eps.eps3 / (AD.k * T3star) - ...
                (i4 + j4) * eps.eps4 / (AD.k * T4star));
            kp = 0;
            if i3 >= 1 && i4 <= i4_max - 1 && j4 <= i4_max - 1
                kp = cached_vv_rate('34_4d', T, i3, i4, j4, AD, opts);
            end
            km = 0;
            if i3 <= i3_max - 1 && i4 >= 1 && j4 >= 1
                kprev = cached_vv_rate('34_4d', T, i3 + 1, i4 - 1, j4 - 1, AD, opts);
                km = detailed_balance_vv(kprev, i3 + 1, i4 - 1, j4 - 1, i3, i4, j4, T, eps, AD, true);
            end
            acc = acc + w * (kp - km);
        end
    end
end
B = acc / den;
end

function [Z3, Z4] = partition_34(T3star, T4star, eps, AD)
T3s = max(T3star, 1e-12);
T4s = max(T4star, 1e-12);
Z3 = 0;
for i3 = 0:(AD.lch4(3) - 1)
    Z3 = Z3 + stat_w3(i3) * exp(-i3 * eps.eps3 / (AD.k * T3s));
end
Z4 = 0;
for i4 = 0:(AD.lch4(4) - 1)
    Z4 = Z4 + stat_w4(i4) * exp(-i4 * eps.eps4 / (AD.k * T4s));
end
end

function k = cached_vv_rate(kind, T, i3, i4, j4, AD, opts)
persistent rate_cache

if ischar(kind) || isstring(kind)
    if strcmpi(string(kind), "clear")
        rate_cache = containers.Map('KeyType', 'char', 'ValueType', 'double');
        k = [];
        return
    end
end

if isempty(rate_cache)
    rate_cache = containers.Map('KeyType', 'char', 'ValueType', 'double');
end

[si, sf, sk, skf, steric, gamma] = vv_forward_states(kind, i3, i4, j4, AD);
key = sprintf('%s_T%.12g_i%d_l%d_j%d_s%.12g_a%.12g_e%.12g_g%.12g_n%d', ...
    kind, T, i3, i4, j4, steric, AD.fho_alpha, AD.fho_e_m, gamma, opts.n_steps_vv);

if opts.use_cache && isKey(rate_cache, key)
    k = rate_cache(key);
    return
end

if abs(gamma - 0.5) < 10 * eps
    k = rates_fho_vv(T, si, sf, sk, skf, steric, AD.fho_alpha, AD.fho_e_m, ...
        'n_steps', opts.n_steps_vv);
else
    k = rates_fho_vv(T, si, sf, sk, skf, steric, AD.fho_alpha, AD.fho_e_m, ...
        'n_steps', opts.n_steps_vv, 'gamma', gamma);
end

if opts.use_cache
    rate_cache(key) = k;
end
end

function [si, sf, sk, skf, steric, gamma] = vv_forward_states(kind, i3, i4, j4, AD)
sk = [0, 0, 0, 0];
skf = [0, 0, 0, 0];
switch lower(kind)
    case '34s'
        si = [0, 0, i3, i4];
        sf = [0, 0, i3 - 1, i4 + 1];
        steric = AD.fho_steric_vv_34s;
        gamma = 0.5;
    case '34d'
        si = [0, 0, i3, i4];
        sf = [0, 0, i3 - 1, i4 + 2];
        steric = AD.fho_steric_vv_34d;
        gamma = AD.fho_gamma_vv34d;
    case '34_4d'
        si = [0, 0, i3, i4];
        sf = [0, 0, i3 - 1, i4 + 1];
        sk = [0, 0, 0, j4];
        skf = [0, 0, 0, j4 + 1];
        steric = AD.fho_steric_vv_34_4d;
        gamma = 0.5;
    otherwise
        error('ch4_relax_terms_simplified:bad_vv_kind', ...
            'Unknown VV process: %s', kind);
end
end

function kb = detailed_balance_vv(kf, i3_i, i4_i, j4_i, i3_f, i4_f, j4_f, T, eps, AD, has_partner)
w_i = stat_w3(i3_i) * stat_w4(i4_i);
w_f = stat_w3(i3_f) * stat_w4(i4_f);
e_i = i3_i * eps.eps3 + i4_i * eps.eps4;
e_f = i3_f * eps.eps3 + i4_f * eps.eps4;
if has_partner
    w_i = w_i * stat_w4(j4_i);
    w_f = w_f * stat_w4(j4_f);
    e_i = e_i + j4_i * eps.eps4;
    e_f = e_f + j4_f * eps.eps4;
end
kb = kf * (w_i / w_f) * exp((e_f - e_i) / (AD.k * T));
end

function sw = stat_w3(level)
sw = 0.5 * (level + 1) * (level + 2);
end

function sw = stat_w4(level)
sw = 0.5 * (level + 1) * (level + 2);
end

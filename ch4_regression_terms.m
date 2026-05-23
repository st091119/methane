function t = ch4_regression_terms(T, varargin)
%CH4_REGRESSION_TERMS  Fast regression-based CH4 vibrational relaxation terms.
%
%   Multi-temperature relaxation source terms evaluated from polynomial
%   regression closures (no FHO rate integration). ~10^4x faster than the
%   state-resolved sums, full-spectrum accurate (closures fit on [9,17,9,20]).
%
%   USAGE
%     t = ch4_regression_terms(T, T13, T24, n)   % 3T -> t.R13, t.R24  [J/(m^3 s)]
%     t = ch4_regression_terms(T, Tv, n)         % 2T -> t.Rvibr (VT only, formulation)
%                                                %       t.Rvibr_hybrid (VT+VV lumped)
%   Inputs:  T   translational-rotational temperature [K]
%            T13,T24 / Tv   vibrational temperature(s) [K]
%            n   number density [m^-3]  (n = p/(k*T))
%   Output struct t also carries per-process R13/R24 in t.proc.<name>.
%
%   WHERE THE MODELS LIVE
%     Coefficients: regression_coefficients.json (9 closures).
%     Reduced fn R = D(driving) * exp(P(polynomial in theta/T)). Reconstructed
%     here with the exact prefactors from regression.tex, e.g.
%        R13^{VVd34} = -n^2 eps3 A34d(T,T13,T24),  R24^{VVd34} = 2 n^2 eps4 A34d.
%     Processes (production params, gamma_vv34d=1.0): VT2, VT4, VV3-4(s),
%     VV3-4(d), VV3-4-4(d), VV3-4-2(d), VV3-2-4(d). The negligible VV3-2(d),
%     VV3-2-2(d) are omitted (sensitivity analysis: >3 orders below dominant).
%
%   EXAMPLE (drop-in fast RHS): see mt_compare_regression.m.

persistent C
if isempty(C), C = load_coeffs(); end

% physical constants / mode quanta [J] (harmonic, from ground state)
h=6.62607015e-34; c=2.99792458e10; kB=1.380649e-23;
omega=[3025.0,1582.7,3156.8,1367.4]; e=omega*h*c;
e2=e(2); e3=e(3); e4=e(4);

if numel(varargin)==3            % --- 3T: (T, T13, T24, n) ---
    T13=varargin{1}; T24=varargin{2}; n=varargin{3};
    [R13v,R24v,names]=rterms(T,T13,T24,n,C,kB,e2,e3,e4);
    t.R13=sum(R13v); t.R24=sum(R24v);
    t.proc=struct(); for k=1:numel(names), t.proc.(names{k})=[R13v(k) R24v(k)]; end
elseif numel(varargin)==2        % --- 2T: (T, Tv, n) ---
    Tv=varargin{1}; n=varargin{2};
    [R13v,R24v,names]=rterms(T,Tv,Tv,n,C,kB,e2,e3,e4);   % all modes at Tv
    t.Rvibr        = R24v(1)+R24v(2);                    % strict 2T: VT2+VT4 only
    t.Rvibr_hybrid = sum(R13v)+sum(R24v);                % hybrid: VT + VV (lumped)
    t.proc=struct(); for k=1:numel(names), t.proc.(names{k})=R13v(k)+R24v(k); end
else
    error('ch4_regression_terms:args','Use (T,T13,T24,n) for 3T or (T,Tv,n) for 2T.');
end
end

% ---- per-process group moments R13, R24 [J/(m^3 s)] (7-process scheme) ----
function [R13v,R24v,names]=rterms(T,T13,T24,n,C,kB,e2,e3,e4)
n2=n^2; u=1000/T; v13=1000/T13; w24=1000/T24; v24=1000/T24; w13=1000/T13; z24=1000/T24;
S2=vtreg(C.vt2,T,T24); S4=vtreg(C.vt4,T,T24);
D34s=1-exp(e3/(kB*T13)-e4/(kB*T24)+(e4-e3)/(kB*T));
D34d=1-exp(e3/(kB*T13)-2*e4/(kB*T24)+(2*e4-e3)/(kB*T));
D324=1-exp(e3/(kB*T13)-e2/(kB*T24)-e4/(kB*T24)+(e2+e4-e3)/(kB*T));
A34s=D34s*exp(poly3(C.a34s,u,v13,w24));
A34d=D34d*exp(poly3(C.a34d,u,v13,w24));
B344=D34d*exp(poly3(C.b344,u,v13,w24));
B342=D324*exp(poly4(C.b342,u,v24,w13,z24));
B324=D324*exp(poly4(C.b324,u,v24,w13,z24));
%       VT2       VT4      VVs34       VVd34         VVd344        VVd342           VVd324
R13v=[0,        0,       -n2*e3*A34s, -n2*e3*A34d,  -n2*e3*B344,  -n2*e3*B342,     -n2*e3*B324];
R24v=[n2*e2*S2, n2*e4*S4, n2*e4*A34s,  2*n2*e4*A34d, 2*n2*e4*B344, n2*(e2+e4)*B342, n2*(e2+e4)*B324];
names={'VT2','VT4','VVs34','VVd34','VVd344','VVd342','VVd324'};
end

function S=vtreg(c,T,Ts)   % VT reduced fn: (1/Ts-1/T)*exp(P(u,v)), 28-coeff deg-6
u=1000/T; v=1000/Ts; P=0; k=1;
for d=0:6, for p=0:d, P=P+c(k)*u^p*v^(d-p); k=k+1; end; end
S=(1/Ts-1/T)*exp(P);
end
function val=poly3(c,u,v,w)   % total-degree-6 polynomial, 84 coeffs
val=0; m=1; for l=0:6, for p=0:l, for q=0:l-p, val=val+c(m)*u^p*v^q*w^(l-p-q); m=m+1; end;end;end
end
function val=poly4(c,u,v,w,z)  % total-degree-6 polynomial, 210 coeffs
val=0; m=1; for l=0:6, for p=0:l, for q=0:l-p, for r=0:l-p-q, val=val+c(m)*u^p*v^q*w^r*z^(l-p-q-r); m=m+1; end;end;end;end
end

function C=load_coeffs()
raw=jsondecode(fileread(fullfile(fileparts(mfilename('fullpath')),'regression_coefficients.json')));
f=fieldnames(raw.coeffs);
for i=1:numel(f), C.(f{i})=raw.coeffs.(f{i})(:); end
end

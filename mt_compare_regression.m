function mt_compare_regression
%MT_COMPARE_REGRESSION  Example: fast 2T vs 3T relaxation using regression models.
%
%   Demonstrates how to run a multi-temperature CH4 relaxation entirely from the
%   regression closures in CH4_REGRESSION_TERMS (no FHO rate integration), on the
%   full vibrational spectrum. Reproduces the thesis figure fig:3t_equal_modes.
%   The whole comparison computes in well under a second.
%
%   The relaxation terms come from one call:
%       t = ch4_regression_terms(T, T13, T24, n);   -> t.R13, t.R24   (3T)
%       t = ch4_regression_terms(T, Tv, n);          -> t.Rvibr        (2T)
%   The ODE state below is the per-molecule mean vibrational energy of each
%   group; temperatures are recovered by inverting <E>(T) and energy conservation
%   3*k*T + <E_vibr> = const.

h=6.62607015e-34; c=2.99792458e10; kB=1.380649e-23;
omega=[3025.0,1582.7,3156.8,1367.4]; eps=omega*h*c; L=[9,17,9,20];

% --- initial conditions for fig:3t_equal_modes (modes start cold/equal, gas hot)
T0=1000; T13_0=300; T24_0=300; Tv0=300; p0=101325;
n0=p0/(kB*T0);

% --- temperature<->energy inversion tables (cheap, no rates)
Tg=linspace(120,1700,500);
E13t=arrayfun(@(T)Ebar(T,[1 3],eps,L,kB),Tg);
E24t=arrayfun(@(T)Ebar(T,[2 4],eps,L,kB),Tg);
Evt =arrayfun(@(T)Ebar(T,[1 2 3 4],eps,L,kB),Tg);
i13=@(E)interp1(E13t,Tg,E,'pchip'); i24=@(E)interp1(E24t,Tg,E,'pchip'); iv=@(E)interp1(Evt,Tg,E,'pchip');

tv=[0, logspace(-12,log10(1),160)]; opts=odeset('RelTol',1e-8,'AbsTol',1e-30);

% --- 3T ---
y0=[Ebar(T13_0,[1 3],eps,L,kB); Ebar(T24_0,[2 4],eps,L,kB)]; C0=3*kB*T0+sum(y0);
[X3,Y3]=ode15s(@(~,y) rhs3(y,C0,n0,i13,i24,kB), tv, y0, opts);
t3=X3(2:end); E13=Y3(2:end,1); E24=Y3(2:end,2);
T13=i13(E13); T24=i24(E24); Tt3=(C0-E13-E24)/(3*kB);

% --- 2T (strict, VT only) and 2T (hybrid, VT+VV lumped) ---
yv0=Ebar(Tv0,[1 2 3 4],eps,L,kB); C0v=3*kB*T0+yv0;
[Xs,Ys]=ode15s(@(~,y) rhs2(y,C0v,n0,iv,kB,false), tv, yv0, opts);
ts=Xs(2:end); Tvs=iv(Ys(2:end)); Tts=(C0v-Ys(2:end))/(3*kB);
[Xh,Yh]=ode15s(@(~,y) rhs2(y,C0v,n0,iv,kB,true), tv, yv0, opts);
th=Xh(2:end); Tvh=iv(Yh(2:end)); Tth=(C0v-Yh(2:end))/(3*kB);

fprintf('Equilibrium T: 3T=%.1f, 2T-strict=%.1f, 2T-hybrid=%.1f K\n', Tt3(end),Tts(end),Tth(end));

f=figure('Visible','off','Position',[100 100 900 560]); hold on;
semilogx(t3,Tt3,'k-','LineWidth',2.2);
semilogx(t3,T13,'-','Color',[0.5 0 0.8],'LineWidth',2.2);
semilogx(t3,T24,'-','Color',[0 0.6 0],'LineWidth',2.2);
semilogx(th,Tth,'r--','LineWidth',1.6);
semilogx(th,Tvh,'--','Color',[0.85 0.4 0],'LineWidth',1.6);
semilogx(ts,Tvs,':','Color',[0.2 0.4 0.9],'LineWidth',1.8);
set(gca,'XScale','log','FontSize',12);
xlabel('t [sec]'); ylabel('Temperature [K]');
legend({'T (3T)','T_{13} (3T)','T_{24} (3T)','T (2T-hybrid)','T_v (2T-hybrid)','T_v (2T-strict, VT only)'},'Location','best');
title(sprintf('2T vs 3T (regression, full spectrum): T_0=%d, T_{13}=T_{24}=%d K',T0,T13_0));
grid on; box on;
if ~exist('imgs','dir'), mkdir('imgs'); end
exportgraphics(f,fullfile('imgs','mt_compare_equal_modes.png'),'Resolution',150);
fprintf('Saved imgs/mt_compare_equal_modes.png\n');
end

function dy=rhs3(y,C0,n0,i13,i24,kB)
T13=max(i13(y(1)),1); T24=max(i24(y(2)),1); T=max((C0-y(1)-y(2))/(3*kB),1);
t=ch4_regression_terms(T,T13,T24,n0);          % <-- regression relaxation terms
dy=[t.R13; t.R24]/n0;
end
function dy=rhs2(y,C0,n0,iv,kB,hybrid)
Tv=max(iv(y),1); T=max((C0-y)/(3*kB),1);
t=ch4_regression_terms(T,Tv,n0);               % <-- regression relaxation terms
if hybrid, dy=t.Rvibr_hybrid/n0; else, dy=t.Rvibr/n0; end
end

function E=Ebar(T,modes,eps,L,kB)   % per-molecule mean energy of the given modes [J]
E=0; for m=modes, E=E+eps(m)*mq(m,T,eps,L,kB); end
end
function mi=mq(mode,T,eps,L,kB)     % mean quantum number <i_m>(T)
Lm=L(mode); em=eps(mode); num=0; den=0;
for i=0:Lm-1, s=sw(i,mode); b=s*exp(-i*em/(kB*T)); den=den+b; num=num+i*b; end
mi=num/den;
end
function s=sw(i,m), switch m, case 1,s=1; case 2,s=i+1; otherwise,s=(i+1)*(i+2)/2; end; end

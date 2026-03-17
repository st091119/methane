function alpha_pop = relative_boltzmann_population_full(T, AD)
% RELATIVE_BOLTZMANN_POPULATION_FULL  Ratio <E_2> / <E_4> of average
%   single-mode vibrational energies at temperature T.
%
%   alpha_pop = relative_boltzmann_population_full(T, AD)

hc    = AD.h * AD.c;
omega = AD.omega;
k_B   = AD.k;

levels2 = AD.lch4(2);
levels4 = AD.lch4(4);

% mode 2
sw2 = zeros(1, levels2);
e2  = zeros(1, levels2);
for n = 0:levels2-1
    sw2(n+1) = stat_weight_mode(n, 2);
    e2(n+1)  = hc * omega(2) * n;
end
b2  = sw2 .* exp(-e2 / (k_B * T));
Z2  = sum(b2);
E2  = sum(e2 .* b2) / Z2;

% mode 4
sw4 = zeros(1, levels4);
e4  = zeros(1, levels4);
for n = 0:levels4-1
    sw4(n+1) = stat_weight_mode(n, 4);
    e4(n+1)  = hc * omega(4) * n;
end
b4  = sw4 .* exp(-e4 / (k_B * T));
Z4  = sum(b4);
E4  = sum(e4 .* b4) / Z4;

alpha_pop = E2 / E4;

end


function sw = stat_weight_mode(level, mode)
switch mode
    case 1, sw = 1;
    case 2, sw = level + 1;
    case 3, sw = 0.5 * (level + 1) * (level + 2);
    case 4, sw = 0.5 * (level + 1) * (level + 2);
    otherwise, error('Invalid mode: %d', mode);
end
end

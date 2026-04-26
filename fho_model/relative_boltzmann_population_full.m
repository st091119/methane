function alpha_pop = relative_boltzmann_population_full(T, AD)
% RELATIVE_BOLTZMANN_POPULATION_FULL  First-excited population ratio for modes 2 and 4.
%
%   alpha_pop = relative_boltzmann_population_full(T, AD)
%
%   alpha = (s_2(1) / s_4(1)) * exp(-(theta_2 - theta_4) / T),
%   where theta_m = h*c*omega_m/k_B.

hc    = AD.h * AD.c;
omega = AD.omega;
k_B   = AD.k;

theta2 = hc * omega(2) / k_B;
theta4 = hc * omega(4) / k_B;
stat_ratio = stat_weight_mode(1, 2) / stat_weight_mode(1, 4);

alpha_pop = stat_ratio .* exp(-(theta2 - theta4) ./ T);

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

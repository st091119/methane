clear

% добавить папку fho_model в путь
addpath(fullfile(fileparts(mfilename('fullpath')), 'fho_model'));

%% входные параметры

% свитчи (используются для одиночного запуска при необходимости)
sw_vr = 'lt';  % правые части: 'lt' или 'sts'
sw_rt = 'm_w'; % ! время релаксации: Landay-Teller: m_w - Milliken-White, vt_rel_time_wang - Wang-Springer, fho - FHO model

% параметры FHO модели (используются при sw_rt = 'fho')
fho_alpha2  = 4.3479e10;     % [м^-1] параметр потенциала Морзе, мода 2
fho_e_m2    = 1434.2;        % [К]    глубина потенциала / k_B, мода 2
fho_alpha4  = 6.6145e10;     % [м^-1] параметр потенциала Морзе, мода 4
fho_e_m4    = 280.0;         % [К]    глубина потенциала / k_B, мода 4
fho_steric2 = 0.0397;        % стерический фактор для моды 2
fho_steric4 = 0.0050;        % стерический фактор для моды 4

% начальные условия
p0 = 101325;        % давление [Па] ! не влияет на решение, будет нужно только для обезразмеривания системы
T0 = 1000;           % температура [К] (какая у нас температура?)
Tv0 = 300;    % колебательная температура [К] для двухтемпературной модели (lt и sts)

% конец интегрирования [с]
t_fin = 1;

% параметры интегрирования
options = odeset('RelTol', 1e-8, 'AbsTol', 1e-8);

%% дополнительные параметры

% загрузка констант
AD = input_data; 

% ! проверка констант
% disp('Input data:');
% disp(AD.klop);
% disp(AD.e1234);
% disp(AD.inds);
disp(AD.stw(1:10));

% ! Время релаксации колебательной энергии CH4-CH4 по формуле Милликена-Уайта
ptau = m_w(T0); % [сек * Па]
fprintf('tau_MW(CH4-CH4) = %.6e сек at %.2f Pa and %.2f K\n', ptau / p0, p0, T0);
ptau_ws = vt_rel_time_wang(T0); % [сек * Па]
fprintf('tau_WS(CH4-CH4) = %.6e сек at %.2f Pa and %.2f K\n', ptau_ws / p0, p0, T0);

% ! Время релаксации по FHO модели
res_fho = relaxation_time_kinetic(T0, fho_alpha2, fho_e_m2, fho_alpha4, fho_e_m4, fho_steric2, fho_steric4);
ptau_fho = res_fho.ptau_total * 101325; % [атм*с] -> [Па*с]
fprintf('tau_FHO(CH4-CH4) = %.6e сек at %.2f Pa and %.2f K\n', ptau_fho / p0, p0, T0);

% ! Распределение Больцмана по колебательным уровням при начальных условиях
% n0 = p0 / (AD.k * T0); % плотность частиц [м^-3]
% ni = vdf(n0, T0, AD);
% disp('Initial vibrational distribution ni:');
% disp(ni);


% среднее время пробега между столкновениями
n0 = p0 / (AD.k * T0); % [м^-3]
sigma0 = pi * AD.r0^2; % [м^2]
tau = (4 * n0 * sigma0 * sqrt(AD.k * T0 / (pi * AD.m)))^(-1); % [сек]
disp(['Mean time between collisions tau = ' num2str(tau) ' sec']);

AD.n0 = n0; AD.T0 = T0; AD.p0 = p0; AD.tau = tau;

%% сохранить свитчи и параметры FHO в структуре AD
AD.sw_rt = sw_rt;

% сохранить параметры FHO в структуре AD (для rpart_mt_lt при sw_rt = 'fho')
AD.fho_alpha2  = fho_alpha2;
AD.fho_e_m2    = fho_e_m2;
AD.fho_alpha4  = fho_alpha4;
AD.fho_e_m4    = fho_e_m4;
AD.fho_steric2 = fho_steric2;
AD.fho_steric4 = fho_steric4;
% для rpart_mt_sts (VV34): при необходимости задайте явно
% AD.fho_steric_vt_v3 = sqrt(fho_steric2 * fho_steric4);
% AD.fho_steric_vv_34s = 0.068;

% интервал интегрирования в безразмерном виде
tspan = [0, t_fin]./tau;

% входной массив начальных условий в безразмерном виде (двухтемпературная модель)
Y0 = [1; Tv0 / T0];
Y0_3t = [1; Tv0 / T0; Tv0 / T0];

%% решение системы для пяти моделей времени релаксации
% LT + MW, LT + Wang-Springer, LT + FHO, STS, трехтемпературная STS
run_cases = struct( ...
    'name',   {'Milliken-White', 'Wang-Springer', 'Landau-Teller (FHO)', 'STS', 'трехтемпературная модель'}, ...
    'rp',     {'rpart_mt_lt',    'rpart_mt_lt',   'rpart_mt_lt',         'rpart_mt_sts', 'rpart_mt_3t_sts'}, ...
    'sw_rt',  {'m_w',            'vt_rel_time_wang', 'fho',              'fho', 'fho'}, ...
    'dim',    {2,                2,               2,                     2,     3} ...
);
results = struct();

for im = 1:length(run_cases)
    AD_run = AD;
    AD_run.sw_rt = run_cases(im).sw_rt;
    RP = str2func(run_cases(im).rp);
    
    fprintf('\n--- Solving with %s ---\n', run_cases(im).name);
    if run_cases(im).dim == 3
        [X_m, Y_m] = ode15s(@(t,y) RP(t, y, AD_run), tspan, Y0_3t, options);
    else
        [X_m, Y_m] = ode15s(@(t,y) RP(t, y, AD_run), tspan, Y0, options);
    end
    
    results(im).time = X_m * tau;
    results(im).T    = Y_m(:,1) * T0;
    if run_cases(im).dim == 3
        results(im).T13 = Y_m(:,2) * T0;
        results(im).T24 = Y_m(:,3) * T0;
    else
        results(im).Tv = Y_m(:,2) * T0;
    end
    results(im).name = run_cases(im).name;
    results(im).dim = run_cases(im).dim;
    
    fprintf('  Final time: %s sec\n', num2str(results(im).time(end)));
    fprintf('  Final T:  %s K\n', num2str(results(im).T(end)));
    if run_cases(im).dim == 3
        fprintf('  Final T13: %s K\n', num2str(results(im).T13(end)));
        fprintf('  Final T24: %s K\n', num2str(results(im).T24(end)));
    else
        fprintf('  Final Tv: %s K\n', num2str(results(im).Tv(end)));
    end
    
    % проверка вычислительной ошибки
    if run_cases(im).dim == 2
        d1 = error_check(Y0 * T0, Y_m(end, :) * T0, AD_run);
        tol = 1e-6;
        if d1 > tol
            fprintf('  Big error: %s\n', num2str(d1));
        else
            fprintf('  Error: %s\n', num2str(d1));
        end
    else
        fprintf('  Error check skipped for 3T model\n');
    end
end

%% графики — все подходы на одном рисунке
colors_T  = {'b', 'r', 'k', [0 0.6 0], [0.5 0 0.8]};
colors_Tv = {'b', 'r', 'k', [0 0.6 0], [0.5 0 0.8]};
styles_T  = {'-', '-', '-', '-', '-'};
styles_Tv = {'--', '--', '--', '--', '--'};

figure; hold on;
legend_entries = {};
for im = 1:length(run_cases)
    semilogx(results(im).time, results(im).T, ...
        [colors_T{im} styles_T{im}], 'LineWidth', 2);
    if results(im).dim == 3
        semilogx(results(im).time, results(im).T13, ...
            [colors_Tv{im} '--'], 'LineWidth', 2);
        semilogx(results(im).time, results(im).T24, ...
            [colors_Tv{im} '-.'], 'LineWidth', 2);
    else
        semilogx(results(im).time, results(im).Tv, ...
            [colors_Tv{im} styles_Tv{im}], 'LineWidth', 2);
    end
    legend_entries{end+1} = ['T — '  results(im).name];
    if results(im).dim == 3
        legend_entries{end+1} = ['T_{13} — ' results(im).name];
        legend_entries{end+1} = ['T_{24} — ' results(im).name];
    else
        legend_entries{end+1} = ['T_v — ' results(im).name];
    end
end
set(gca, 'XScale', 'log');
xlabel('t [sec]');
ylabel('Temperature [K]');
legend(legend_entries, 'Location', 'best');
title(['T_0 = ' num2str(T0) ' K, T_{v0} = ' num2str(Tv0) ' K']);
grid on;

%% ══════════════════════════════════════════════════════════════
%  Сравнение трёх моделей: Milliken-White, Wang-Springer, FHO
%  ══════════════════════════════════════════════════════════════
T_comp = linspace(200, 1300, 80);
ptau_comp_mw  = zeros(size(T_comp));
ptau_comp_ws  = zeros(size(T_comp));
ptau_comp_fho = zeros(size(T_comp));

for iT = 1:length(T_comp)
    Ti = T_comp(iT);
    ptau_comp_mw(iT)  = m_w(Ti);               % [Па*с]
    ptau_comp_ws(iT)  = vt_rel_time_wang(Ti);   % [Па*с]
    res_i = relaxation_time_kinetic(Ti, fho_alpha2, fho_e_m2, fho_alpha4, fho_e_m4, fho_steric2, fho_steric4);
    ptau_comp_fho(iT) = res_i.ptau_total * 101325; % [атм*с] -> [Па*с]
end

figure;
semilogy(T_comp, ptau_comp_mw / 101325, 'b-', 'LineWidth', 2); hold on;
semilogy(T_comp, ptau_comp_ws / 101325, 'r--', 'LineWidth', 2);
semilogy(T_comp, ptau_comp_fho / 101325, 'k-.', 'LineWidth', 2);
xlabel('T [K]');
ylabel('p\tau [atm \cdot s]');
legend('Milliken-White', 'Wang-Springer', 'FHO');
title('Сравнение времён VT релаксации CH4-CH4');
grid on;



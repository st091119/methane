function dy = rpart_mt_lt(~, y, AD)
% функция расчета правых частей на основе соотношения Ландау-Теллера

%% переменные в безразмерном виде
T_b = y(1);
Tv_b = y(2);

% число неизвестных
num_var = length(y);

%% переменные в размерном виде
% температура [К]
temp = T_b * AD.T0;
tv = Tv_b * AD.T0;    

%% расчет колебательных распределений (векторизовано на упрощенном наборе индексов)

% inds = AD.inds;    % [N x 4]
% I0 = inds(:,1); J0 = inds(:,2); K0 = inds(:,3); L0 = inds(:,4);

% Веса и энергии
w = AD.stw;        % [N x 1]
E = AD.e1234;      % [N x 1]

% Фактор
% sumE = I0 .* AD.e1000 + J0 .* AD.e0100 + K0 .* AD.e0010 + L0 .* AD.e0001; % [Дж]
% fac = -(sumE) ./ (AD.k * tv);                                            % [безр.]
fac = -(E) ./ (AD.k * tv);                                             % [безр.]

% Статсумма и колебательная энергия газа при температуре tv
z_vibr = sum(w .* exp(fac));
% e_vibr = (1 / (z_vibr * AD.m)) * sum(w .* E .* fac); % на единицу массы [Дж/кг]

e_sum_e_sum =  sum(w .* E ./ (AD.k * tv) .* exp(fac))^2; % безр.
e_sum_square = sum(w .* (E ./ (AD.k * tv)).^2 .* exp(fac)); % безр.

% e_sum_e_sum =  sum(w .* sumE ./ (AD.k * AD.T0) .* exp(fac)) * sum(w .* sumE ./ (AD.k * tv) .* exp(fac)); % безр.
% e_sum_square = sum(w .* sumE .* sumE ./ (AD.k * AD.T0) / (AD.k * tv) .* exp(fac)); % безр.

%% релаксационные члены
% обратная величина времени релаксации [сек^-1]
if strcmp(AD.sw_rt, 'fho')
    res_fho = relaxation_time_kinetic(temp, AD.fho_alpha2, AD.fho_e_m2, ...
                                      AD.fho_alpha4, AD.fho_e_m4, ...
                                      AD.fho_steric2, AD.fho_steric4);
    p_tau = res_fho.ptau_total * 101325; % [атм*с] -> [Па*с]
else
    fun_times = str2func(AD.sw_rt);
    p_tau = fun_times(temp);
end

times_inv = 1 / (p_tau / AD.p0);

% функции расчета rho*E_m/n [Дж]
mEv = @(t) sum(w .* E .* exp(-E ./ (AD.k * t))) / sum(w .* exp(-E ./ (AD.k * t)));
% disp(['mEv(temp) = ' num2str(mEv(temp)) ' J, temp = ' num2str(temp) ' K']);
% disp(['mEv(tv) = ' num2str(mEv(tv)) ' J, tv = ' num2str(tv) ' K']);
% mEv = @(t) sum(w .* sumE .* exp(-sumE ./ (AD.k * t))) / sum(w .* exp(-sumE ./ (AD.k * t)));

% размерные релаксационные члены [Дж/сек]
RVIBR_DIM = (mEv(temp) - mEv(tv)) * times_inv; 
% disp(['RVIBR_DIM = ' num2str(RVIBR_DIM) ' J/s']);

% безразмерные релаксационные члены
kT0 = AD.k * AD.T0;
RVIBR = RVIBR_DIM * AD.tau / kT0;

%% составляем матрицу коэффициентов перед производными А
% единичная матрица 
A = eye(num_var);

% уравнение сохранения энергии
% T
A(1, 1) = 3;

% Tv
A(1, 2) = e_sum_square / (z_vibr) - e_sum_e_sum / (z_vibr^2);

% уравнение сохранения кол. энергии
A(2, 2) = A(1, 2);

AA = sparse(A);

%% составляем вектор-столбец правых частей B
B = zeros(num_var,1);
B(2) = RVIBR;

dy = AA^(-1) * B;

end
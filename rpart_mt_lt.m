function dy = rpart_mt_lt(~, y, AD)
% функция расчета правых частей на основе соотношения Ландау-Теллера


%% переменные в безразмерном виде
var = num2cell(y);
[T_b, Tv_b] = deal(var{:});

% число неизвестных
num_var = length(y);

%% переменные в размерном виде
% температура [К]
temp = T_b*AD.T0;
tv = Tv_b*AD.T0;    

%% дополнительные параметры
kT0 = AD.k*AD.T0; % [Дж]

% колебательная энергия молекул в безр. виде
e1234_b = AD.e1234 ./ kT0;

% вспомогательные величины
xi = -(AD.e1234)./ (AD.k * tv);

% статистические сумма
Zv = sum(AD.stw .* exp(xi));

% суммы в уравнениях:
si_exi = AD.stw .* exp(xi);
S_ei_si_xi_exi = sum((AD.e1234 + AD.e0000) .* xi .* si_exi);

S_ei_si_exi = sum((AD.e1234 + AD.e0000) .* si_exi);

S_si_xi_exi = sum(xi .* si_exi);


%% релаксационные члены
% обратная величина времени релаксации [сек^-1] по Милликену-Уайту
times_inv = m_w(AD.p0, temp);

% функции расчета rho*E_m/n [Дж]
mE = @(t) sum(AD.stw .* AD.e1234 .* exp(-AD.e1234/(AD.k*t))) ./ sum(AD.stw .* exp(-AD.e1234./(AD.k*t)));

% размерные релаксационные члены [Дж/сек]
RVIB = (mE(temp) - mE(tv)) * times_inv; 

% безразмерные релаксационные члены
RVIB = RVIB * AD.tau/kT0;

%% составляем матрицу коэффициентов перед производными А
% единичная матрица 
A = eye(num_var);

% уравнение сохранения энергии
% T
A(1,1) = 2.5;

% Tv
A(1,2) = -S_ei_si_xi_exi/(Zv*Tv_b) + S_ei_si_exi*S_si_xi_exi/(Zv^2*Tv_b);


% уравнение сохранения кол. энергии
A(2,2) = -S_ei_si_xi_exi/(Zv*Tv_b) + S_ei_si_exi*S_si_xi_exi/(Zv^2*Tv_b);

AA = sparse(A);

%% составляем вектор-столбец правых частей B
B = zeros(num_var,1);
B(2) = RVIB;

dy = AA^(-1)*B;

end
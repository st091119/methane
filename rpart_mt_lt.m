function dy = rpart_mt_lt(~, y, AD)
% функция расчета правых частей на основе соотношения Ландау-Теллера

%% переменные в безразмерном виде
T_b = y(1);
Tv_b = y(2);

% число неизвестных
num_var = length(y);

%% переменные в размерном виде
% температура [К]
temp = T_b*AD.T0;
tv = Tv_b*AD.T0;    
%% расчет колебательных распределений
[ii,jj,kk, ll] = ind2sub(AD.lch4, AD.ind_1d_in_3d);
i_state = ii-1; j_state = jj-1; k_state = kk-1; l_state = ll-1;

E2 = i_state*AD.e1234(2, 1, 1, 1) / AD.k;
E3 = j_state*AD.e0100 / AD.k;
E4 = k_state*AD.e0010 / AD.k;
E5 = l_state*AD.e0001 / AD.k;

AD.lch4 = [9, 17, 9, 20];
AD.d = [1, 2, 3, 3]; %степени вырожденности по модам
z_vibr = 0;
e_vibr = 0;
for l = 1:AD.lch4(4)
    for k = 1:AD.lch4(3)
        for j = 1:AD.lch4(2)
            for i = 1:AD.lch4(1) 
                z_vibr = z_vibr + AD.stw(i, j, k, l) * (-(i * AD.e1000 + j * AD.e0100 + k * AD.e0010 + l * AD.e0001)/(AD.k * tv));
                e_vibr = 0;
            end
        end
    end
end


expE = exp(-E2-E3-E4-E5);
Zv = sum(AD.stw .* expE);

%% дополнительные параметры
kT0 = AD.k*AD.T0; % [Дж]

% колебательная энергия молекул в безр. виде
eco2i_b = AD.eco2_i/kT0;
eco20_b = AD.eco2_0/kT0;

e12_b = AD.e12/kT0;
e3_b = AD.e3/kT0;

% вспомогательные величины
% xi в уравнения будут отличаться в зависимости от суммирования
% уравнение внутренней энергии:
xi12_IE = -(E2+E3)/t12;
xi3_IE = -E4/t3;
xi = xi12_IE + xi3_IE;

% уравнения колебательных энергий:
xi12_VE = -AD.e12/(AD.k*t12);
xi3_VE = -AD.e3/(AD.k*t3);

% статистические суммы объединенной и антисимметричной мод
Z12 = sum(AD.stw12 .* exp(xi12_VE));
Z3 = sum(exp(xi3_VE));

% суммы в уравнении внутренней энергии:
si_exi = AD.stw .* exp(xi);

S_ei_si_xi12_exi = sum((eco2i_b+eco20_b) .* xi12_IE .* si_exi);
S_ei_si_xi3_exi = sum((eco2i_b+eco20_b) .* xi3_IE .* si_exi);

S_ei_si_exi = sum((eco2i_b+eco20_b) .* si_exi);

S_si_xi12_exi = sum(xi12_IE .* si_exi);
S_si_xi3_exi = sum(xi3_IE .* si_exi);

% суммы в уравнениях колебательных энергий:
si_e12_exi12 = AD.stw12 .* e12_b .* exp(xi12_VE);
e3_exi3 = e3_b .* exp(xi3_VE);

S_si_e12_xi12_exi12 = sum(si_e12_exi12 .* xi12_VE);
S_e3_xi3_exi3 = sum(e3_exi3 .* xi3_VE);

S_si_e12_exi12 = sum(si_e12_exi12);
S_e3_exi3 = sum(e3_exi3);

S_si_xi12_exi12 = sum(AD.stw12 .* xi12_VE .* exp(xi12_VE));
S_xi3_exi3 = sum(xi3_VE .* exp(xi3_VE));

%% релаксационные члены
% обратная величина времени релаксации [сек^-1] по Милликену-Уайту
times_inv = m_w(AD.p0, temp);

% функции расчета rho*E_m/n [Дж]
mE12 = @(t) sum(AD.stw12 .* AD.e12 .* exp(-AD.e12/(AD.k*t))) / sum(AD.stw12 .* exp(-AD.e12/(AD.k*t)));
mE3 = @(t) sum(AD.e3 .* exp(-AD.e3/(AD.k*t))) / sum(exp(-AD.e3/(AD.k*t)));

% размерные релаксационные члены [Дж/сек]
RVIB = (mEv(temp) - mEv(tv)) * times_inv; 

% безразмерные релаксационные члены
RVIB = RVIB * AD.tau/kT0;

%% составляем матрицу коэффициентов перед производными А
% единичная матрица 
A = eye(num_var);

% уравнение сохранения энергии
% T
A(1,1) = 2.5;

% Tv12
A(1,2) = -S_ei_si_xi12_exi/(Zv*T12_b) + S_ei_si_exi*S_si_xi12_exi/(Zv^2*T12_b);

% Tv3
A(1,3) = -S_ei_si_xi3_exi/(Zv*T3_b) + S_ei_si_exi*S_si_xi3_exi/(Zv^2*T3_b);

% уравнение сохранения кол. энергии
% объединенная мода
A(2,2) = -S_si_e12_xi12_exi12/(Z12*T12_b) + S_si_e12_exi12*S_si_xi12_exi12/(Z12^2*T12_b);

% антисимметричнаяя мода
A(3,3) = -S_e3_xi3_exi3/(Z3*T3_b) + S_e3_exi3*S_xi3_exi3/(Z3^2*T3_b);

AA = sparse(A);

%% составляем вектор-столбец правых частей B
B = zeros(num_var,1);
B(2) = RVIB_12;
B(3) = RVIB_3;

dy = AA^(-1)*B;

end
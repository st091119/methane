function calculate_vibrational_energies()
    % Константы
    h = 6.62607015e-34;     % Постоянная Планка, Дж·с
    c = 2.99792458e10;      % Скорость света, см/с
    
    % Данные из таблицы (моды колебаний молекулы метана)
    modes = [
        1, 3025.0;  % Mode 1: Symmetrical
        2, 1582.7;  % Mode 2: Twisting  
        3, 3156.8;  % Mode 3: Antisymmetrical
        3, 1367.4   % Mode 4: Scissoring
    ];
    
    % Вычисление энергий для различных квантовых состояний
    fprintf('Энергии колебательных состояний молекулы метана:\n\n');
    
  % вектор колебательной энергии молекул CH4 в расчете от 0-го уровня [Дж]
    e0001 = calculate_epsilon(h, c, modes, 0, 0, 0, 1);
    e0010 = calculate_epsilon(h, c, modes, 0, 0, 1, 0);
    e0100 = calculate_epsilon(h, c, modes, 0, 1, 0, 0);
    e1000 = calculate_epsilon(h, c, modes, 1, 0, 0, 0);
    e0000 = calculate_epsilon(h, c, modes, 0, 0, 0, 0);

    fprintf('ε₀₀₀₁ =%.4f  Дж\n', e0001);
    
    
    fprintf('ε₀₀₁₀ = %.4f Дж\n', e0010);
    
    
    fprintf('ε₀₁₀₀ = %.4f Дж\n', e0100);
    
   
    fprintf('ε₁₀₀₀ = %.4f Дж\n', e1000);
    
   
    fprintf('ε₀₀₀₀ = %.4f Дж\n', e0000);
end

function epsilon = calculate_epsilon(h, c, modes, i1, i2, i3, i4)
    % ε = hc * Σ[ω_m * (i_m + d_m/2)]
    
    i_m = [i1, i2, i3, i4];
    
    % Вычисление суммы
    sum_val = 0;
    for m = 1:4
        d_m = modes(m, 1);      % Степень вырождения
        omega_m_e = modes(m, 2); % Спектроскопическая постоянная, см⁻¹
        sum_val = sum_val + omega_m_e * (i_m(m) + d_m/2);
        print(sum_val);
    end
     epsilon = h * c * sum_val;
    end
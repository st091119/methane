function d1 = error_check(Y0, Y1, AD)
    % Y0: [T0, Tv0]
    % Y1: [T_vec, Tv_vec] - матрицы из ODE
    
    T0 = Y0(1); TV0 = Y0(2);
    
    % Начальная энергия
    ni0 = vdf(1, TV0, AD);
    E0 = (1.5 + 1.5)*AD.k*T0 + sum(AD.e1234 .* ni0);

    num_steps = size(Y1, 1);
% невязка
    u1 = zeros(num_steps, 1);
    
    for i = 1:num_steps
        Ti = Y1(i, 1);
        TVi = Y1(i, 2);
        
        ni = vdf(1, TVi, AD);
        Ei = (1.5 + 1.5)*AD.k*Ti + sum(AD.e1234 .* ni);
        u1(i) = Ei - E0;
    end

    d1 = max(abs(u1) / E0);
end
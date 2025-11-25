function dat = train

%% константы
dat.c = 29979245800;         % скорость света [см/сек]
dat.h = 6.62607015e-34;    % постоянная Планка [Дж*сек]
dat.hbar = dat.h/(2*pi);
dat.Na = 6.02214076e23;    % постоянная Авогадро [1/моль]
dat.R = 8.314462618;       % универсальная газовая постоянная [Дж/моль/К]
dat.k = 1.380649e-23;      % постоянная Больцмана [Дж/К]
dat.Torr = 133.322;        % [Торр] -> [Па]

% колебательная частота молекулы CH4 по модам (см^(-1))
omega = [3025.0, 1582.7, 3156.8, 1367.4];

dat.lch4 = [9, 17, 9, 20];
dat.d = [1, 2, 3, 3]; %степени вырожденности по модам
stw = zeros (dat.lch4(1), dat.lch4(2), dat.lch4(3), dat.lch4(4));
e1234 = zeros (dat.lch4(1), dat.lch4(2), dat.lch4(3), dat.lch4(4));
for l = 1:dat.lch4(4)
    for k = 1:dat.lch4(3)
        for j = 1:dat.lch4(2)
            for i = 1:dat.lch4(1) 
                ishka = [i, j, k, l];
                stw(i,j, k, l) = (j * k * (k + 1) * l * (l + 1)) / 4;
                e1234(i, j, k, l) = 0;
                for m = 1:4
                    e1234(i,j, k, l) = e1234(i,j,k,l) + (dat.h * dat.c * ((ishka(m) - 1) *(stw(m) + dat.d(m) / 2)));
                end
            end
        end
    end
end
fprintf('e1234(1,1,1,1) = %.6e Дж\n', e1234(1,1,1,1));
fprintf('e1234(2,1,1,1) = %.6e Дж\n', e1234(2,1,1,1));
fprintf('e1234(2,2,2,2) = %.4e Дж\n', e1234(2,2,2,2));

fprintf('stw(1,1,1,1) = %.4f\n', stw(1,1,1,1));
fprintf('stw(2,1,1,1) = %.4f\n', stw(2,1,1,1));
fprintf('stw(2,2,2,2) = %.4f\n', stw(2,2,2,2));

dat.stw = stw;
dat.e1234 = e1234;
dat.omega = omega;
end

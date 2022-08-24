Mcoarse = load('B:\TDOAmodel_200m.mat');
MfineFolder = 'B:\modelFiles_10mFrom200m';

load('D:\SOCAL_E_63\xwavTables\instrumentLocs.mat')  % calculated in D:\MATLAB_addons\gitHub\wheresWhaledo\experiments\calcSigma.m
h = [0,0,0; h];

N = 100;

for n = 1:N
    [errCoarse(n), Icoarse] = min(sqrt(sum((Mcoarse.wloc - ).^2, 2)));
    TDOAcoarse(n, :) = Mcoarse.TDOA(Icoarse, :);

    modelFile = ['TDOAmodel_10m_n=', num2str(Icoarse, '%05.f')];
    Mfine = load(fullfile(MfineFolder, modelFile));
    [errFine(n), Ifine] = min(sqrt(sum((Mfine.wloc - w(n, :)).^2, 2)));
    TDOAfine(n, :) = Mfine.TDOA(Ifine, :);
end

figure(1)
for n = 1:6
    subplot(6,1,n)
    plot(TDOAcoarse(:, n), '.')
    hold on
    plot(TDOAfine(:, n), '.')
    hold off
end

figure(2)
for n = 1:6
    subplot(6,1,n)
    plot(TDOAcoarse(:, n+6), '.')
    hold on
    plot(TDOAfine(:, n+6), '.')
    hold off
end

figure(3)
for n = 1:6
    subplot(6,1,n)
    plot(TDOAcoarse(:, n+12), '.')
    hold on
    plot(TDOAfine(:, n+12), '.')
    hold off
end

figure(4)
scatter3(h(:, 1), h(:, 2), h(:, 3), 'rs')
hold on
scatter3(w(:, 1), w(:, 2), w(:, 3), 'filled')
hold off
% MAIN  Script di sperimentazione per il progetto: nodi di Leja approssimati
clear; clc; close all;

% Mesh per i nodi (scegli 1e4 o 1e5)
Mmesh = 1e5;                      
xmesh = linspace(-1, 1, Mmesh).'; 
deg_max = 50;                     
xeval = linspace(-1, 1, 5000).';  
f = @(x) 1./(x - 1.3);            

if deg_max + 1 > numel(xmesh)
    error('La mesh ha %d punti ma serve almeno deg_max+1 = %d.', numel(xmesh), deg_max+1);
end

times1 = zeros(deg_max,1);    
times2 = zeros(deg_max,1);    
Leb1   = zeros(deg_max,1);    
Leb2   = zeros(deg_max,1);    
err_leja = zeros(deg_max,1);
err_equi = zeros(deg_max,1);

for d = 1:deg_max
    tic; z1 = DLP(xmesh, d); times1(d) = toc;
    tic; z2 = DLP2(xmesh, d); times2(d) = toc;

    Leb1(d) = leb_con(z1, xeval);
    Leb2(d) = leb_con(z2, xeval);

    z_leja = z2(:);                       
    z_equi = linspace(-1, 1, d+1).';
    f_leja = f(z_leja);
    f_equi = f(z_equi);

    j = 0:d;
    V_leja = cos(acos(z_leja) * j);
    V_equi = cos(acos(z_equi) * j);

    c_leja = V_leja \ f_leja;   
    c_equi = V_equi \ f_equi;

    V_eval = cos(acos(xeval) * j);
    p_leja = V_eval * c_leja;
    p_equi = V_eval * c_equi;

    err_leja(d) = max(abs(p_leja - f(xeval)));
    err_equi(d) = max(abs(p_equi - f(xeval)));

    fprintf('d=%2d | t1=%.4fs, t2=%.4fs | Leb1=%.3e Leb2=%.3e | errLeja=%.3e errEqui=%.3e\\n', ...
        d, times1(d), times2(d), Leb1(d), Leb2(d), err_leja(d), err_equi(d));
end

% ---- Grafici ----
scriptDir = fileparts(mfilename('fullpath'));
imgDir = fullfile(scriptDir, '..', 'doc', 'img');
if ~exist(imgDir,'dir'), mkdir(imgDir); end

figure;
plot(1:deg_max, times1, 'o-', 'DisplayName', 'DLP (produttoria)'); hold on;
plot(1:deg_max, times2, 's-', 'DisplayName', 'DLP2 (LU Chebyshev)');
xlabel('Grado d'); ylabel('Tempo [s]'); grid on; legend('Location','northwest');
title(sprintf('Tempi computazionali (Mmesh=%d)', Mmesh));
exportgraphics(gcf, fullfile(imgDir,'tempi.png'), 'Resolution', 300);

figure;
semilogy(1:deg_max, Leb2, 's-', 'DisplayName', 'Leja (DLP2)'); hold on;
semilogy(1:deg_max, Leb1, 'o-', 'DisplayName', 'Leja (DLP)');
grid on;
xlabel('Grado d'); ylabel('Costante di Lebesgue (semilog)');
title('Costante di Lebesgue per nodi di Leja');
legend('Location','northwest');
exportgraphics(gcf, fullfile(imgDir,'lebesgue.png'), 'Resolution', 300);

figure;
semilogy(1:deg_max, err_leja, 's-', 'DisplayName', 'Leja (DLP2)'); hold on;
semilogy(1:deg_max, err_equi, 'o-', 'DisplayName', 'Equispaziati'); grid on;
xlabel('Grado d'); ylabel('Errore massimo su [-1,1]');
title('Confronto accuratezza interpolante (base di Chebyshev)');
legend('Location','southwest');
exportgraphics(gcf, fullfile(imgDir,'errori.png'), 'Resolution', 300);

fprintf('Figure salvate in %s: tempi.png, lebesgue.png, errori.png\n', imgDir);

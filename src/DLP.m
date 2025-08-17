function dlp = DLP(x, d)
% DLP  Discrete Leja Points (Algoritmo 1: massimizzazione della produttoria)
%   dlp = DLP(x, d)
%   INPUT:
%     x : vettore (colonna o riga) con i punti della mesh in [-1,1]
%     d : grado del polinomio -> produce d+1 nodi
%   OUTPUT:
%     dlp : vettore riga con i d+1 punti di Leja approssimati
%
%   NOTE:
%   - Il primo nodo e' x(1).
%   - Ogni nodo successivo massimizza la produttoria dei moduli delle distanze
%     dai nodi gia' scelti, valutata sui soli punti della mesh x.
%
x = x(:);                % Colonna
M = numel(x);
if d >= M
    error('DLP: d deve essere < length(x).');
end

% Primo nodo: x(1)
dlp = zeros(1, d+1);
dlp(1) = x(1);

available = true(M,1);
available(1) = false;
prodvals = ones(M,1);

for s = 2:d+1
    dist = abs(x - dlp(s-1));
    dist(~available) = 1;
    prodvals = prodvals .* dist;

    prodvals(~available) = -inf;
    [~, idx] = max(prodvals);
    dlp(s) = x(idx);
    available(idx) = false;
end
end

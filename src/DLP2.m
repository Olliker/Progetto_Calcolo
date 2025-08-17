function dlp = DLP2(x, d)
% DLP2  Discrete Leja Points (Algoritmo 2: LU pivoting su Vandermonde di Chebyshev)
%   dlp = DLP2(x, d)
%   INPUT:
%     x : vettore (colonna o riga) con i punti della mesh in [-1,1]
%     d : grado del polinomio -> produce d+1 nodi
%   OUTPUT:
%     dlp : vettore riga con i d+1 punti di Leja approssimati
%
%   Costruisce la matrice di tipo Vandermonde con base di Chebyshev:
%     V_{i,j} = cos((j-1)*arccos(x_i))
%   Poi esegue LU con pivoting parziale.
%   In piu', forziamo che il PRIMO nodo sia x(1), come richiesto.
%
x = x(:);                
M = numel(x);
if d >= M
    error('DLP2: d deve essere < length(x).');
end

xt = min(max(x, -1), 1);
theta = acos(xt);
j = 0:d;
V = cos(theta * j);

[~, ~, p] = lu(V, 'vector');   

% Forza x(1) come primo nodo
if p(1) ~= 1
    p = [1; p(p~=1)];
end

idx = p(1:d+1);
dlp = x(idx).';
end

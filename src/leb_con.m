function L = leb_con(z, x)
% LEB_CON  Approssimazione della costante di Lebesgue sui punti x
%   L = leb_con(z, x)
%   INPUT:
%     z : vettore riga dei nodi dell'interpolante (length n = d+1)
%     x : vettore colonna dei punti in [-1,1] dove valutare la funzione di Lebesgue
%   OUTPUT:
%     L : scalare reale, max_{x} lambda_n(x)

z = z(:).';           
x = x(:);             
n = numel(z);

% Pesi baricentrici
W = ones(1, n);
for i = 1:n
    diffs = z(i) - z;
    diffs(i) = 1;
    W(i) = 1 / prod(diffs);
end

X = x * ones(1, n);
Z = ones(numel(x), 1) * z;
A = W ./ (X - Z);

tol = 1e-14;
is_node = abs(X - Z) < tol;

S = sum(A, 2);
lambda = sum(abs(A), 2) ./ abs(S);

any_node = any(is_node, 2);
lambda(any_node) = 1;

L = max(lambda);
end

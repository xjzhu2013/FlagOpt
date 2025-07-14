function [f,G] = nleigflag(Y,dim,M)
[n, nd] = size(Y);
d = length(dim);
dim = [0, dim];
G = zeros(n, nd);
f = 0;
for i = 1:d
    Yi = Y(:, dim(i)+1:dim(i+1));
    MYi = M*Yi;
    tri = sum(sum(Yi.*MYi));
    f = f - tri^2;
    G(:, dim(i)+1:dim(i+1)) = -4*tri*MYi;
end
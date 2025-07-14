function P = Proj(D,Y,dim)
[n, nd] = size(Y);
d = length(dim);
dim = [0, dim];
P = zeros(n, nd);
for i = 1:d
    Yi = Y(:, dim(i)+1:dim(i+1));
    Di = D(:, dim(i)+1:dim(i+1));
    YitDi = Yi'*Di;
    P(:, dim(i)+1:dim(i+1)) = Di - (Yi*(YitDi-YitDi') + Y*(D'*Yi));
end
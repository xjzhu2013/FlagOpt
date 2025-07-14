function [Q,B] = preExp(Z,Y)
[n, nd] = size(Y);
[U, ~] = eigs(eye(n) - Y*Y', n - nd);
Q = [Y, U];
B1 = Q'*Z;
B2 = B1(nd+1:n,:);
B = [B1, [-B2'; zeros(n-nd)]];
B = (B - B')/2;
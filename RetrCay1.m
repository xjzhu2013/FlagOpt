function [Yt,M] = RetrCay1(t,U,VtU,VtY,Y)
nd = size(Y,2);
M = linsolve(eye(2*nd) - 0.5*t*VtU, VtY);
Yt = Y + t*U*M;
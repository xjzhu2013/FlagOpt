function [Yt,M5,M6] = RetrCay2(t,Z,YtZ,ZtZ,YYtZ,ZtYYtZ,M1,M2,M3,M4,Y)
nd = size(Y,2);
M5 = Y + 0.5*t*M1;
M6 = YtZ' + 0.5*t*M2;
A = eye(nd) + 0.5*t*M6;
B = t*ZtZ - 0.5*t*ZtYYtZ + 0.25*t^2*M4 + YtZ';
Yt = Y + t*Z - 0.5*t*YYtZ + 0.25*t^2*M3 - 0.5*t*M5*linsolve(A,B); 
function [YtZ,ZtZ,YYtZ,ZtYYtZ,M1,M2,M3,M4] = preCay2(Z,Y)
YtZ = Y'*Z;
ZtZ = Z'*Z;
YYtZ = Y*YtZ;
ZtYYtZ = YtZ'*YtZ;
M1 = Z - YYtZ;
M2 = ZtZ - ZtYYtZ;
M3 = M1*YtZ;
M4 = M2*YtZ;
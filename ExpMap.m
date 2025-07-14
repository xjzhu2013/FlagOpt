function [Yt,Qt] = ExpMap(t,Q,B,nd)
Qt = Q*expm(t*B);
Yt = Qt(:,1:nd);
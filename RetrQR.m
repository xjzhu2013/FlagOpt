function Yt = RetrQR(t,Z,Y)
[Yt, R] = qr(Y + t*Z, 0);
Yt = Yt * diag(sign(sign(diag(R))+.5));
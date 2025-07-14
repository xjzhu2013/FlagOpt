function T = VecTran1(t,U,VtU,VtY,M)
T = U*(VtY + t*VtU*M);
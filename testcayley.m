function outavrg = testcayley(Y, opts, dim, M)
outavrg.iter = 0;
outavrg.fval = 0;
outavrg.nrmg = 0;
outavrg.feasi = 0;
for k = 1:10
    Ys = Y(:,(k-1)*100+1:k*100);
    [~, out] = FlagCGcay(@nleigflag, Ys, opts, dim, M);
    outavrg.iter = outavrg.iter + out.iter;
    outavrg.fval = outavrg.fval + out.fval;
    outavrg.nrmg = outavrg.nrmg + out.nrmg;
    outavrg.feasi = outavrg.feasi + out.feasi;
end
outavrg.iter = outavrg.iter/10;
outavrg.fval = outavrg.fval/10;
outavrg.nrmg = outavrg.nrmg/10;
outavrg.feasi = outavrg.feasi/10;
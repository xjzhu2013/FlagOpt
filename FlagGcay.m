function [Y, out] = FlagGcay(fun, Y0, opts, dim, varargin)
% Gradient method on Flag manifold: Cayley transform
% Input: fun - objective 
%        Y0 - initial point
%        opts - options
%        dim - dimensions of flag
%        varargin - data in objective
% Output: Y - solution
%         out - output information

% problem options
c = opts.c;
t = opts.t;
rho = opts.rho;
gtol = opts.gtol;
mxitr = opts.mxitr;
methodretr = opts.methodretr;

% initialization
nd = size(Y0,2);
Y = Y0;

% compute function value and gradient at Y
[f, G] = feval(fun, Y, dim, varargin{:});
Grad = Proj(G, Y, dim);
nrmg = norm(Grad,'fro');

% record of norm of gradient
recg = nrmg;

% stopping tolerance
tol = gtol*nrmg;

% main loop
for k = 1:mxitr
    
    % convergence test
    if nrmg <= tol
        break
    end
    
    % pre-Cayley transform
    if methodretr == 1
        [U,VtU,VtY] = preCay1(-Grad,Y);
    else
        [YtZ,ZtZ,YYtZ,ZtYYtZ,M1,M2,M3,M4] = preCay2(-Grad,Y);
    end

    % backtracking line search
    count = 0;
    while count < 5
        if methodretr == 1
            [Yt,~] = RetrCay1(t,U,VtU,VtY,Y);
        else
            [Yt,~,~] = RetrCay2(t,-Grad,YtZ,ZtZ,YYtZ,ZtYYtZ,M1,M2,M3,M4,Y);
        end
        [ft, Gt] = feval(fun, Yt, dim, varargin{:});
        ratio = (f - ft)/(t*nrmg^2);
        if ratio >= c
            break
        else
            count = count + 1;
            t = t*rho;
        end
    end
    
    %re-orthogonalization for Yt
    YtY = Yt'*Yt;
    feasi = norm(YtY - eye(nd), 'fro'); 
    if feasi > 1e-13
        Yt = Yt*(YtY)^(-1/2);
        %[U, ~, V] = svd(Yt, 0); Yt = U*V'; 
        [ft, Gt] = feval(fun, Yt, dim, varargin{:});
    end
    
    % new gradient and its norm
    Gradt = Proj(Gt, Yt, dim);
    nrmg = norm(Gradt,'fro');
    
    % BB step
    S = t*Grad;
    W = Gradt - Grad;
    %W = Gt - G;
    if mod(k,2) == 0
        t = sum(sum(S.*S))/abs(sum(sum(S.*W)));
    else
        t = abs(sum(sum(S.*W)))/sum(sum(W.*W));
    end
    %t = sum(sum(S.*S))/abs(sum(sum(S.*W)));
    t = max(min(t, 1e20), 1e-20);
    
    % update Y, function value and gradient at Y
    Y = Yt;
    f = ft;
    %G = Gt;
    Grad = Gradt;
    
    % update record
    recg_new = [recg; nrmg];
    recg = recg_new;
    
end

% output
out.iter = k;
out.fval = f;
out.nrmg = nrmg;
out.recg = recg;
out.feasi = feasi;
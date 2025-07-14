function [Y, out] = FlagGexp(fun, Y0, opts, dim, varargin)
% Gradient method on Flag manifold: Exponential map
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
    
    % pre-exponential map
    [Q,B] = preExp(-Grad,Y);

    % backtracking line search
    count = 0;
    while count < 5
        [Yt,~] = ExpMap(t,Q,B,nd);
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
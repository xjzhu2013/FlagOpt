function [Y, out] = FlagCGdec(fun, Y0, opts, dim, varargin)
% Conjugate gradient method on Flag manifold: Polar or QR decomposition
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
methodcg = opts.methodcg;

% initialization
nd = size(Y0,2);
Y = Y0;

% compute function value and gradient at Y
[f, G] = feval(fun, Y, dim, varargin{:});
Grad = Proj(G, Y, dim);
nrmg = norm(Grad,'fro');

% search direction and directional derivative
Z = -Grad;
df = -nrmg^2;

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

    % line search
    count = 0;
    while count < 5
        if methodretr == 1
            Yt = RetrPol(t,Z,Y);
        else
            Yt = RetrQR(t,Z,Y);
        end
        [ft, Gt] = feval(fun, Yt, dim, varargin{:});
        ratio = (f - ft)/(-t*df);
        if ratio >= c
            break
        else
            count = count + 1;
            
            % backtracking
            t = t*rho;
            
            % interpolation
            %t = -0.5*df*t^2/(ft - f - df*t);
            
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
    
    % new gradient with its norm
    Gradt = Proj(Gt, Yt, dim);
    nrmgt = norm(Gradt,'fro');
    
    % difference of gradient
    W = Gradt - Grad;
    
    % vector transport
    T = Proj(Z, Yt, dim);
    
    % update beta
    switch methodcg
        case 'FR' % Fletcher-Reeves
            beta = (nrmgt/nrmg)^2;
        case 'PR' % Polak-Ribiere+
            beta = sum(sum(Gradt.*W))/nrmg^2;
            beta = max(beta, 0);
        case 'HS' % Hestenes-Stiefel+
            beta = sum(sum(Gradt.*W))/sum(sum(Z.*W));
            beta = max(beta, 0);
        case 'DY' % Dai-Yuan+
            beta = nrmgt^2/sum(sum(Z.*W));
            beta = max(beta, 0);
    end
    
    % decay beta
    beta = beta/k;
    
    % new search direction and directional derivative
    Zt = -Gradt + beta*T;    
    df = sum(sum(Gradt.*Zt)); 
    if df > 0
        Zt = -Gradt;
        df = -nrmgt^2;
    end
    
    % BB step
    if mod(k,2) == 0
        t = t*sum(sum(Z.*Z))/abs(sum(sum(Z.*W)));
    else
        t = t*abs(sum(sum(Z.*W)))/sum(sum(W.*W));
    end
    %t = t*sum(sum(Z.*Z))/abs(sum(sum(Z.*W)));
    t = max(min(t, 1e20), 1e-20); 
    
    % update Y, function value and gradient at Y
    Y = Yt;
    f = ft;
    Grad = Gradt;
    nrmg = nrmgt;
    Z = Zt;
    
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
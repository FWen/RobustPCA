function [L,S,out] = lq_lq_l2_bcd(Y,lamda,q1,q2,para);
% lq_lq_l2_bcd solves
%
%   minimize ||L||_Sq1^Sq1 + lamda*||S||_q2^q2 + \beta/2*||Y-L-S||_F^2
%
% Inputs 
%	1=>q1,q2=>0
%	Ltrue, Strue: for debug, for calculation of errors
%   L0,S0: intialization
% Outputs
%	L,S: the recovery
%	out.el, out.es: the error with respect to the true

[m, n] = size(Y);

%Initialize
if(nargin==5 & isfield(para, 'MAX_ITER'))
    MAX_ITER = para.MAX_ITER;
else
    MAX_ITER = 2e3;
end

if(nargin==5 & isfield(para, 'TOL'))
    ABSTOL = para.TOL;
else
    ABSTOL = 1e-6;
end

if(nargin==5 & isfield(para, 'L0'))
    L = para.L0;
else
    L = zeros(m,n);
end

if(nargin==5 & isfield(para, 'S0'))
    S = para.S0;
else
    S = zeros(m,n);
end
    
mu  = norm(Y); 
phi = 1e-5;

out.el = [];

for iter = 1:MAX_ITER

    Lm1 = L; 
    Sm1 = S;
	
    % for acceleration of the algorithm
    mu = mu * 0.97;
    
    % L-update 
    [U,V,D] = svd( (Y-S+mu*phi*L) / (1+mu*phi) );
    v = shrinkage_Lq(diag(V), q1, 1, 1/mu+phi);
    indx = find(v);
    L = U(:,indx)*diag(v(indx))*D(:,indx)';
    
    % S-update
    T = (Y-L+mu*phi*S)/(1+mu*phi);
    S = reshape(shrinkage_Lq(T(:), q2, lamda, (1/mu+phi)), m, n);  
            
    %Check for convergence
    if (norm(L-Lm1,'fro')< sqrt(m*n)*ABSTOL) & (norm(S-Sm1,'fro')< sqrt(m*n)*ABSTOL)
        num_stop = num_stop + 1;
        if num_stop==3
            break;
        end
    else
        num_stop = 0;
    end
    
end

end

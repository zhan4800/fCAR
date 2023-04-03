function Rivmat=autocor(autoc)
%function [Rmat,Rivmat, Lmat,R1,L1,Rinv]=autocor(autoc)
% This function transforms the auto correlation parameters to the autocorrelation matrix R.
% Input: autoc a structure that contains the following:
%        autoc.rho: correlation parameters.
%        autoc.type: the type of the autocorrelation. 'AR1',
%        'Mixed.p','DampedExponential'
%        autoc.Ti: the time points on which each measurements are
%                          obtained.
% Output:
% R1: the cell of correlation marices of all measurement sets.
% L1: the cell of cholesky decomposition of R1.
% Rmat: the diagonal block matrix formed by R1 matices.
% Lmat: the diagonal block matrix formed by L1 matices.
Nt=length(autoc.Ti); % the number of blocks in the R matrix.
R1=cell(Nt,1);
L1=cell(Nt,1);
Rinv=cell(Nt,1);
rho=autoc.rho;

for i = 1:Nt
    ti = autoc.Ti{i};
    nti = length(ti);     
    switch autoc.type
        case 'AR1' % In this case, rho is a scalar.
           R1{i} = rho.^abs(repmat(ti,1,nti)-repmat(ti',nti,1)); % ti is assumed to be a column.           
           %L1{i} = chol(R1{i},'lower');
        case 'Mixed.p'
           %ri = rho{i};
           pi = length(rho);  % In case of Mixed.p, rho is a cell of length Nt. In each cell, length(rho{i}) is the order.
           R1{i} = eye(nti);
           for j = 1: pi
              R1{i} = R1{i}+rho(j)*diag(ones(nti-j,1),j)+rho(j)*diag(ones(nti-j,1),-j); 
           end              
           %L1{i} = chol(R1{i},'lower');
        case 'DampedExponential' % In case of Damped exponential, I need rho, ts, where corr(X(t1),X(t2)=rho^{|t1-t2|^ts}).          
           ts = autoc.ts; % ts is a scalar in [0,1]
           R1{i}=rho.^(abs(repmat(ti,1,nti)-repmat(ti',nti,1)).^ts);                      
    end
    
    if any(eig(R1{i})<=0) 
        Rinv{i} = pinv(R1{i});
    else
        L1{i} = chol(R1{i},'lower');
        Linv=solve_tril(L1{i},eye(nti));
        Rinv{i}=Linv'*Linv;
    end
end
%Rmat=blkdiag(R1{:});
%Lmat=blkdiag(L1{:});
Rivmat=blkdiag(Rinv{:});
        

% 
% switch autoc.type
%     
%     case 'AR1'
%         autoc.rho
%         
%     case 'Mixed1'
%         err=laprnd(n,wavespecs.K,zeros(n,wavespecs.K), repmat(sqrt(S0)',n,1));
%         u=laprnd(M,wavespecs.K,zeros(M,wavespecs.K),repmat(sqrt(Q0)',M,1));
%     case 'Mixedp'
%         nu=3;
%         err=trnd(nu,n,wavespecs.K).*repmat(sqrt(S0)'/sqrt(nu/(nu-2)),n,1);
%         u=trnd(nu,M,wavespecs.K).*repmat(sqrt(Q0)'/sqrt(nu/(nu-2)),M,1);
%     case 'DampExpo'
%         nu=2;
%         err=trnd(nu,n,wavespecs.K).*repmat(sqrt(S0)'/2,n,1);
%         u=trnd(nu,M,wavespecs.K).*repmat(sqrt(Q0)'/2,M,1);    
% end



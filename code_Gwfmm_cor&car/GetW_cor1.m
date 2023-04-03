function W=GetW_cor1(model,D,Rinv)
% This function get the W matrix for the case that correlation of E_{jk}^*
% is R.
%%%%% W=GetW(model,D)
%%%%%   Compute cross-products matrices
%%%%%           X'RinvX       X'RinvZ       X'RinvD
%%%%%           Z'RinvX       Z'RinvZ       Z'RinvD
%%%%%           D'RinvX       D'RinvZ       D'RinvD
% Input:    
%
%           model = structure with integer elements:
%               X = (n x p) design matrix for fixed effects functions;
%               Z = cell array containg H matrices (each n x m_h); design
%                       matrices for random effects functions; 
%           D = (n x K) matrix of wavelet coefficients
%
% Output:   W = structure with elements:
%           XtX,XtZ,XtD,ZtD,DtD
W.XtX=model.X'*Rinv*model.X;
W.XtD=model.X'*Rinv*D;
W.DtD=D'*Rinv*D;

H=model.H;
if (H>0)
    Z=model.Z{1};
    if (H>1)
        for h=2:H
            Z=[Z,model.Z{h}];    %#ok<*AGROW> %% concatenate columns of model.Z{h}, h=1,...,H to form single Z matrix.
        end
    end
else
    Z=0;       %#zhu#% may be make Z to be zeros(n,1) will make matrix product straight forward. 
end
W.ZtZ=Z'*Rinv*Z;

if Z==0
    W.ZtD=zeros(1,size(D,2)); 
    W.XtZ=zeros(size(model.X,2),1);  
else
    W.ZtD=Z'*Rinv*D;
    W.XtZ=model.X'*Rinv*Z;
end


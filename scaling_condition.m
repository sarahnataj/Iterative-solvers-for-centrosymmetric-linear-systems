% Centrosymmetry-preserving row and column equilibration algorithm
% Given nonsingular n\times n matrix A, this program computes
% nonsingular diagonal matrices R and S such that B = RAS has the property that
% max_k | b_ik|  = max_k | b_ki|  = 1 for all i and R = S if A = A^T . 
%
% Input: The 1D Poisson and biharmonic spectral differentiatian matrices,
%       tol is a convergence tolerance.
% Output: The graphs of 1D Poisson and biharmonic matrices after equilibration.
% Author: Sarah Nataj, email:sarah.nataj@gmail.com
clear all
close all
nn=5:1:30;
tol=1e-15;
for ii=1:size(nn,2)
    n=nn(ii);
    %A=PseudoSpectral1D(n+1,2);
    A=PseudoSpectral1D(n+1,4);
    n=size(full(A),1);nn(ii)=n;
    cond1(ii)=cond(A);
    R=eye(n);S=eye(n);r=zeros(n,1);s=zeros(n,1);
    while (max(abs(r-1))>tol)||(max(abs(s-1))>tol)
        for i=1:n
            r(i)=1/sqrt(norm(A(i,:),inf));
            s(i)=1/sqrt(norm(A(:,i),inf));
        end
        A=diag(r)*A*diag(s);
        R=diag(r)*R;
        S=S*diag(s);
    end
    B=A;
    cond2(ii)=cond(B);
end
loglog(nn,cond1,'-o',nn,cond2,'-*',nn,0.25*nn.^4,'--k','LineWidth',2)
legend('condition(A)', 'condition(RAS)','O(N^4)','Location','NorthWest')
ylabel('Condition number');
title('1D biharmonic equation')
xlabel('N')
plotformat(1.5,6)
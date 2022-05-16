function [R,S]=scaling(A)
% Centrosymmetry-preserving row and column equilibration algorithm
% Given nonsingular n\times n matrix A, this algorithm computes
% nonsingular diagonal matrices R and S such that B = RAS has the property that
% max_k | b_ik|  = max_k | b_ki|  = 1 for all i and R = S if A = A^T . 
% (Algorithm 3)
% Input:A is a given square matrix
%       tol is a convergence tolerance.
% Output: Diagonal matrices R and S such that B=RAS is the equilibrated matrix
%        with smaller condition number with respect to A.
% Author: Sarah Nataj, email:sarah.nataj@gmail.com
n=size(A,1);
tol=1e-15;
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

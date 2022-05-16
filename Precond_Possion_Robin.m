function Precond_Possion_Robin(n,precase,tol,toltp)
% PCG preconditioned by incomplete double-cone factorization
% for 2D Poisson operator with Robin boundary conditions
% Input: n is the number of the collocation nodes
%        precase: two preconditioners:
%               precase=1: the preconditioner is the centrosymmetric part 
%               of the original matrix
%               precase=2 the approximated preconditioner from solving the
%               same PDE with Robin boundary conditions with average a.
% 
%        toltp = ilutp drop tolerance,
%        tol = GMRES tolerance.
% Output: The number of iterations for GMRES and density of proposed preconditioners.
% Author: Sarah Nataj, email:sarah.nataj@gmail.com
%
[A,M]=PoissonRobin2D(n);% A is orginal matrix, M is approximted one with average a.
%r=eig(full(A));%min(r), max(r)
x0 = zeros(size(A,1),1);
b=randn(size(A,1),1);
maxit = size(A,1);
if precase==1
    disp('Preconditioner is centrosymetric part of the matrix')
    J=flipud(speye(size(A,1)));
    M=0.5*(A+J*A*J);
else
    disp('Preconditioner is the approximated matrix from solving the same PDE with average a.')
end
%s=eig(full(M\A));%min(s), max(s)
disp('Condition number of the matrix+++++++++++++++++++++++++++++');
condition=cond(full(A))
L0=chol(M,'lower');
density_chol=(nnz(L0)+nnz(L0')-size(M,1))/nnz(M)
disp('PCG without preconditioner+++++++++++++++++++++++++++++++++');
[x1,fl1,rr1,it1,rv1] = pcg(A,b,tol,maxit,[],[],x0);
it_non=it1
%
disp('PCG with ICHOL(0) factor as preconditioner++++++++++++++++++++')
L=ichol(M,struct('shape','lower'));
[x2,fl2,rr2,it2,rv2] = pcg(A,b,tol,maxit,L,L',x0);
density_ichol=(nnz(L)+nnz(L')-size(M,1))/nnz(M)
it_ichol=it2
%
disp('PCG with IXX(0) factor as preconditioner++++++++++++++++++++++')
X=ixx(M);
[x3,fl3,rr3,it3,rv3] = pcg(A,b,tol,maxit,X,X',x0);
density_ixx=(nnz(X)+nnz(X')-size(M,1))/nnz(M)
it_ixx=it3
%
disp('PCG with ICHOLT fator as preconditioner++++++++++++++++++++')
Lt=ichol(M,struct('shape','lower','type','ict','droptol',toltp,'michol','on'));
density_ict=(nnz(Lt)+nnz(Lt')-size(M,1))/nnz(M)
[x4,fl4,rr4,it4,rv4] = pcg(A,b,tol,maxit,Lt,Lt',x0);
iter_icholt=it4

%
disp('PCG with IXXT factor as preconditioner++++++++++++++++++++++')
Xt=ixxt(M,toltp);
density_ixt=(nnz(Xt)+nnz(Xt')-size(M,1))/nnz(M)
[x5,fl5,rr5,it5,rv5] = pcg(A,b,tol,maxit,Xt,Xt',x0);
iter_ixxt=it5

figure(1); semilogy(0:length(rv2)-1,rv2/norm(b),'-s',0:length(rv3)-1,rv3/norm(b),'-*',...
    0:length(rv4)-1,rv4/norm(b),'->b', 0:length(rv5)-1,rv5/norm(b),'-<r',...
    'LineWidth',1.15)
yline(tol,'r--')
legend('ICHOL(0)','IXX(0)','ICHOLT(10^{-3})','IXXT(10^{-3})','tolerance','Location','NorthEast')
xlabel('Iteration number')
ylabel('Relative residual')
title('Robin boundary conditions');

% figure(2); plot(r,0,'*r'), title('unpreconditioned')
%
% figure(3); plot(s,0,'*b'), title('preconditioned')

plotformat(1.5,6)






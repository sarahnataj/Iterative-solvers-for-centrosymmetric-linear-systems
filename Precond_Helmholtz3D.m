function Precond_Helmholtz3D(n,a,tol,toltp)
% GMRES preconditioned by incomplete double-cone factorization
% for 3D  Helmholtz equation -(\Delta +a)u=f
% Here, A is  spectral differentiation matrix associated with Helmholtz equation.
% Input: n is the number of the collocation nodes
%       toltp = ilutp drop tolerance
%       tol = GMRES tolerance
%       a is square of the wave number
% Output: the number of iterations for GMRES and density of proposed preconditioners.
% Author: Sarah Nataj, email:sarah.nataj@gmail.com

[A, uex]=Helmholtz3D(n+1,a);
condition=cond(full(A))
x0 = zeros(size(A,1),1);
b=randn(size(A,1),1);
maxit = size(A,1);
[L, U]=lu(A);
density_lu=(nnz(L)+nnz(U)-size(A,1))/nnz(A)
disp('GMRES without preconditioner+++++++++++++++++++++++++++++++');
[x0,fl0,rr0,it0,rv0] = gmres(A,b,[],tol,maxit,[],[],x0);
iter_non=it0(2)
disp('GMRES with ILU(0) factors as preconditioner+++++++++++++++++++');
[L1, U1]=ilu(A);
density_ilu=(nnz(L1)+nnz(U1)-size(A,1))/nnz(A)
[x1,fl1,rr1,it1,rv1] = gmres(A,b,[],tol,maxit,L1,U1);
iter_ilu=it1(2)
%
disp('GMRES with IXY(0) factors as preconditioner++++++++++++++++++++')
[X1, Y1]=ixy(A);
density_ixy=(nnz(X1)+nnz(Y1)-size(A,1))/nnz(A)
[x2,fl2,rr2,it2,rv2]= gmres(A,b,[],tol,maxit,X1,Y1);
iter_ixy=it2(2)
disp('GMRES with ILUTP factors as preconditioner+++++++++++++++++');
[L2, U2,P]=ilu(A,struct('type','ilutp','droptol',toltp));
density_ilutp=(nnz(L2)+nnz(U2)-size(A,1))/nnz(P*A)
[x3,fl3,rr3,it3,rv3]= gmres(P*A,P*b,[],tol,maxit,L2,U2);
iter_ilutp=it3(2)
disp('GMRES with IXYTP factors as preconditioner++++++++++++++++++')
[X2, Y2, P1]=ixytp(A,toltp);
density_ixytp=(nnz(X2)+nnz(Y2)-size(A,1))/nnz(P1*A)
[x4,fl4,rr4,it4,rv4]= gmres(P1*A,P1*b,[],tol,maxit,X2,Y2);
iter_ixytp=it4(2)


figure (1); semilogy(0:length(rv1)-1,rv1/norm(b),'-s',...
    0:length(rv2)-1,rv2/norm(b),'-*',...
    0:length(rv3)-1,rv3/norm(b),'->b',...
    0:length(rv4)-1,rv4/norm(b),'-<r',...
    'LineWidth',1.15)
yline(tol,'r--')
plotformat(1.5,6)
legend('ILU(0)','IXY(0)','ILUTP(10^{-3})','IXYTP(10^{-3})','tolerance','Location','NorthEast')
xlabel('Iteration number')
ylabel('Relative residual')
title('3D Helmholtz equation with Dirichlet boundary conditions');
plotformat(1.5,6)

end


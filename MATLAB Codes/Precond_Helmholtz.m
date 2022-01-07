% GMRES preconditioned by incomplete double-cone factorization
% for 3D  Helmholtz equation -(\delta +a)u=f
% Here, A is 3D second order spectral differentiation matrix associated
% with Helmholtz equation.
% Input: n is the number of the collocation nodes
%       toltp = ilutp drop tolerance
%       tol = GMRES tolerance
%       a is square of the wave number (here we choose close to eigenvalues
%       of the 3D discrete Laplacian spectral operator with multiplicity 6 
%       for the case n=20)
% Output: the number of iterations for GMRES and density of proposed preconditioners.
% Author: Sarah Nataj, email:sarah.nataj@gmail.com
clc
clear all;
close all;
n=20;
tol = 1e-6;
toltp=[1e-4];
%a=[219.5988, 605.2396,1052.7138, 2195.9064, 2434.9594,3362.1876, ...
%4096.6722, 9794.2739, 18759.0058];
a=[219.5988];
for i=1:size(a,2)
    A=PseudoSpectral_Helmholtz3D(n+1,a(i));
    condition(i)=cond(full(A));
    x0 = zeros(size(A,1),1);
    b=randn(size(A,1),1);
    maxit = size(A,1);
    %disp('without preconditioner+++++++++++++++++++');
    [L, U]=lu(A);
    density_lu(i)=(nnz(L)+nnz(U)-size(A,1))/nnz(A);
    [x0,fl0,rr0,it0,rv0] = gmres(A,b,[],tol,maxit,[],[],x0);
    %fl0
    iter(i)=it0(2);
    %disp('ilu+++++++++++++++++++');
    [L1, U1]=ilu(A);
    density_ilu(i)=(nnz(L1)+nnz(U1)-size(A,1))/nnz(A);
    [x1,fl1,rr1,it1,rv1] = gmres(A,b,[],tol,maxit,L1,U1);
    %fl1
    iter_ilu(i)=it1(2);
    %
    %disp('ixy++++++++++++++++++++++++++++++++++++++');
    [X1, Y1]=ixy(A);
    density_ixy(i)=(nnz(X1)+nnz(Y1)-size(A,1))/nnz(A);
    [x2,fl2,rr2,it2,rv2]= gmres(A,b,[],tol,maxit,X1,Y1);
    %fl2
    iter_ixy(i)=it2(2);
    %disp('ilutp+++++++++++++++++++++++++++++++++++');
    for ii=1:size(toltp,2)
        [L2, U2,P]=ilu(A,struct('type','ilutp','droptol',toltp(ii)));
        density_ilutp(i,ii)=(nnz(L2)+nnz(U2)-size(A,1))/nnz(P*A);
        [x3,fl3,rr3,it3,rv3]= gmres(P*A,P*b,[],tol,maxit,L2,U2);
        iter_ilutp(i)=it3(2);
    end
    %disp('ixytp+++++++++++++++++++++++++++++++++++');
    for ii=1:size(toltp,2)
        [X2, Y2, P1]=ixytp(A,toltp(ii));
        density_ixytp(i,ii)=(nnz(X2)+nnz(Y2)-size(A,1))/nnz(P1*A);
        [x4,fl4,rr4,it4,rv4]= gmres(P1*A,P1*b,[],tol,maxit,X2,Y2);
        iter_ixytp(i)=it4(2);
    end
end

plot(a,iter_ilu,'-s',a,iter_ixy,'-*',...
    a,iter_ilutp,'->b',a,iter_ixytp,'-<r','LineWidth',1.15);
legend('ILU(0)','IXY(0)','ILUTP(10^{-3})','IXYTP{10^{-3}}','Location','NorthWest')
title('3D Helmholtz equation');
xlabel('\lambda^2')
ylabel('iterations')

disp('Square of considered wave numbers:')
a
disp('Condition number of associated matrix')
condition
disp('Iteration number for GMRES without preconditioner')
iter
disp('Iteration number for GMRES with ilu factors as preconditioner')
iter_ilu
disp('Iteration number for GMRES with ixy factors as preconditioner')
iter_ixy
disp('Iteration number for GMRES with ilutp factors as preconditioner')
iter_ilutp
disp('Iteration number for GMRES with ixytp factors as preconditioner')
iter_ixytp
disp('Density of factors')
density_ilu
density_ixy
density_ilutp=density_ilutp';density_ilutp,
density_ixytp=density_ixytp';density_ixytp,




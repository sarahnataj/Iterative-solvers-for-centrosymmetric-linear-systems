% GMRES preconditioned by incomplete double-cone factorization
% for 2D/3D Poisson equation and
%     2D Poisson equation with variable coefficients (diffusion equation).
% Here, A is 2D/3D second order spectral differentiation matrices.
% Input: n is the number of the collocation nodes
%       toltp = ilutp drop tolerance
%       tol = GMRES tolerance
% Output: The number of iterations for GMRES and density of proposed preconditioners.
% Author: Sarah Nataj, email:sarah.nataj@gmail.com
clc
clear all;
close all;
n=11;
tol = 1e-6;
toltp=[1e-3,1e-4];
%A=PseudoSpectral2D(n,2);pcase=1;disp('2D Poisson operator');
A=PseudoSpectral3D(n);pcase=2;disp('3D Poisson operator');
%A=PseudoSpectral_Poisson_varcoef(n,2,10);pcase=3;disp('2D Poisson operator with variable cofficients')
x0 = zeros(size(A,1),1);
b=randn(size(A,1),1);
maxit = size(A,1);
disp('Condition number of the matrix+++++++++++++++++++++++++++++');
condition=condest(A)
disp('GMRES without preconditioner+++++++++++++++++++++++++++++++');
[L, U]=lu(A);
density_lu=(nnz(L)+nnz(U)-size(A,1))/nnz(A)
[x0,fl0,rr0,it0,rv0] = gmres(A,b,[],tol,maxit,[],[],x0);
it0(2),
%
disp('GMRES with ilu factors as preconditioner+++++++++++++++++++');
[L1, U1]=ilu(A);
density_ilu=(nnz(L1)+nnz(U1)-size(A,1))/nnz(A)
[x1,fl1,rr1,it1,rv1] = gmres(A,b,[],tol,maxit,L1,U1);
it1(2),
%
disp('GMRES with ixy factors as preconditioner++++++++++++++++++++')
[X1, Y1]=ixy(A);
density_ixy=(nnz(X1)+nnz(Y1)-size(A,1))/nnz(A)
[x2,fl2,rr2,it2,rv2]= gmres(A,b,[],tol,maxit,X1,Y1);
it2(2),
%
disp('GMRES with ilutp factors as preconditioner+++++++++++++++++');
for ii=1:size(toltp,2)
[L2, U2,P]=ilu(A,struct('type','ilutp','droptol',toltp(ii)));
density_ilutp(ii)=(nnz(L2)+nnz(U2)-size(A,1))/nnz(P*A);
[x3,fl3,rr3,it3,rv3]= gmres(P*A,P*b,[],tol,maxit,L2,U2);
iter_ilutp(ii)=it3(2);
end
density_ilutp
iter_ilutp
%
disp('GMRES with ixytp factors as preconditioner++++++++++++++++++')
for ii=1:size(toltp,2)
[X2, Y2, P1]=ixytp(A,toltp(ii));
density_ixytp(ii)=(nnz(X2)+nnz(Y2)-size(A,1))/nnz(P1*A);
[x4,fl4,rr4,it4,rv4]= gmres(P1*A,P1*b,[],tol,maxit,X2,Y2);
iter_ixytp(ii)=it4(2);
end
density_ixytp
iter_ixytp
%
semilogy(0:length(rv1)-1,rv1/norm(U1\(L1\(P*b))),'-s',...
   0:length(rv2)-1,rv2/norm(Y1\(X1\b)),'-*',...
   0:length(rv3)-1,rv3/norm(U2\(L2\(P1*b))),'->b',...
   0:length(rv4)-1,rv4/norm(Y2\(X2\b)),'-<r', 'LineWidth',1.15);
yline(tol,'r--');
legend('ilu(0)','ixy(0)','ilutp','ixytp','Tolerance','Location','NorthEast')
if  pcase==1
    title('2D Poisson equation');
   elseif pcase==2
    title('3D Poisson equation');
elseif pcase==3
title('2D variable-cofficient Poisson equation')
end
xlabel('Iteration number')
ylabel('Relative residual')
plotformat(1.5,6)



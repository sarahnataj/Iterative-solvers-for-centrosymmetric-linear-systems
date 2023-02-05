function Precond_Poisson(n,ecase,tol,toltp)
% GMRES preconditioned by incomplete double-cone factorization
% for 2D/3D Poisson equation and
%     2D Poisson equation with variable coefficients (diffusion equation).
% Here, A is associated spectral differentiation matrices.
% Input: n = the number of the collocation nodes
%       ecase= the equation
%       toltp = ilutp drop tolerance
%       tol = GMRES tolerance
% Output: The number of iterations for GMRES and density of proposed preconditioners.
% Author: Sarah Nataj, email:sarah.nataj@gmail.com
%
collocation_nodes=n+1
if  ecase==1
    A=PseudoSpectral2D(n,2);
   elseif ecase==2
    A=PseudoSpectral3D(n);
elseif ecase==3
   A=Poisson_varcoef(n,2,10);
end
x0 = zeros(size(A,1),1);
b=randn(size(A,1),1);
maxit = size(A,1);
size_A=size(A)
disp('Condition number of the matrix+++++++++++++++++++++++++++++');
condition_A=cond(full(A))
[L, U]=lu(A);
density_lu=(nnz(L)+nnz(U)-size(A,1))/nnz(A)
disp('GMRES without preconditioner+++++++++++++++++++++++++++++++');
[x0,fl0,rr0,it0,rv0] = gmres(A,b,[],tol,maxit,[],[],x0);
it_non=it0(2),
%
disp('GMRES with ILU(0) factors as preconditioner+++++++++++++++++++');
[L1, U1]=ilu(A);
density_ilu=(nnz(L1)+nnz(U1)-size(A,1))/nnz(A)
[x1,fl1,rr1,it1,rv1] = gmres(A,b,[],tol,maxit,L1,U1);
it_ilu=it1(2),
%
disp('GMRES with IXY(0) factors as preconditioner++++++++++++++++++++')
[X1, Y1]=ixy(A);
density_ixy=(nnz(X1)+nnz(Y1)-size(A,1))/nnz(A)
[x2,fl2,rr2,it2,rv2]= gmres(A,b,[],tol,maxit,X1,Y1);
it_ixy=it2(2),
%
disp('GMRES with ILUTP factors as preconditioner+++++++++++++++++');
[L2, U2,P]=ilu(A,struct('type','ilutp','droptol',toltp));
density_ilutp=(nnz(L2)+nnz(U2)-size(A,1))/nnz(P*A)
[x3,fl3,rr3,it3,rv3]= gmres(P*A,P*b,[],tol,maxit,L2,U2);
it_ilutp=it3(2)
%
disp('GMRES with IXYTP factors as preconditioner++++++++++++++++++')
[X2, Y2, P1]=ixytp(A,toltp);
density_ixytp=(nnz(X2)+nnz(Y2)-size(A,1))/nnz(P1*A)
[x4,fl4,rr4,it4,rv4]= gmres(P1*A,P1*b,[],tol,maxit,X2,Y2);
it_ixytp=it4(2)
%
semilogy(0:length(rv1)-1,rv1/norm(U1\(L1\(P*b))),'-s',...
   0:length(rv2)-1,rv2/norm(Y1\(X1\b)),'-*',...
   0:length(rv3)-1,rv3/norm(U2\(L2\(P1*b))),'->b',...
   0:length(rv4)-1,rv4/norm(Y2\(X2\b)),'-<r', 'LineWidth',1.15);
yline(tol,'r--');
legend('ILU(0)','IXY(0)', 'ILUTP(10^{-4})', 'IXYTP(10^{-4})','tolerance','Location','NorthEast')
if  ecase==1
    title('2D Poisson equation');
   elseif ecase==2
    title('3D Poisson equation');
elseif ecase==3
title('2D diffusion equation')
end
xlabel('Iteration number')
ylabel('Relative residual')
plotformat(1.5,6)

end

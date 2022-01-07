% PCG preconditioned by incomplete double-cone factorization
% for 2D/3D Poisson operator
% Here A is 2D second order SPD spectral differentiation matrix using
% Legendre Lobatto points. We also assume SPD matrices arising from solving
% 2D/3D Poisson equation using finite differences schemes.
%
% Input: n is the number of collocation nodes,
%       toltp = ilutp drop tolerance,
%       tol = GMRES tolerance.
%Output: The number of iterations for GMRES and density of proposed preconditioners.
clc;
clear all;
close all;
tol = 1e-6;
toltp=[1e-3,1e-4];
n=11;pcase=1;%symmetric Legendre spectral matrix for 2D Poisson
%n=30;pcase=2;%finite difference for 2D Poisson
%n=10;pcase=3;%finite difference for 3D Poisson
if pcase==1
    A=SymLeg(n+1);
elseif pcase==2
    h=1/n;A=toeplitz([2; -1; zeros(n-3,1)])/(h^2);
    A=kron(speye(n-1),A)+kron(A,speye(n-1));
elseif pcase==3
    h=1/n;A=toeplitz([2; -1; zeros(n-3,1)])/(h^2);
    A=kron(kron(speye(n-1),A),speye(n-1))+...
        kron(kron(A,speye(n-1)),speye(n-1))...
        +kron(kron(speye(n-1),speye(n-1)),A);
end
x0 = zeros(size(A,1),1);
b=randn(size(A,1),1);
maxit = size(A,1);
disp('Condition number of the matrix+++++++++++++++++++++++++++++');
condition=condest(A)
disp('PCG without preconditioner+++++++++++++++++++++++++++++++++');
L0=chol(A,'lower');
density_chol=(nnz(L0)+nnz(L0')-size(A,1))/nnz(A)
[x1,fl1,rr1,it1,rv1] = pcg(A,b,tol,maxit,[],[],x0);
it1
%
disp('PCG with ichol factor as preconditioner++++++++++++++++++++')
L=ichol(A,struct('shape','lower'));
[x2,fl2,rr2,it2,rv2] = pcg(A,b,tol,maxit,L,L',x0);
density_ichol=(nnz(L)+nnz(L')-size(A,1))/nnz(A)
it2
%
disp('PCG with ixx factor as preconditioner++++++++++++++++++++++')
X=ixx(A);
[x3,fl3,rr3,it3,rv3] = pcg(A,b,tol,maxit,X,X',x0);
density_ixx=(nnz(X)+nnz(X')-size(A,1))/nnz(A)
it3
%
disp('PCG with icholt fator as preconditioner++++++++++++++++++++')
for ii=1:size(toltp,2)
    Lt=ichol(A,struct('shape','lower','type','ict','droptol',toltp(ii),'michol','off'));
    density_ict(ii)=(nnz(Lt)+nnz(Lt')-size(A,1))/nnz(A);
    [x4,fl4,rr4,it4,rv4] = pcg(A,b,tol,maxit,Lt,Lt',x0);
    iter_ict(ii)=it4;
end
density_ict
iter_ict
%
disp('PCG with ixxt factor as preconditioner++++++++++++++++++++++')
for ii=1:size(toltp,2)
    Xt=ixxt(A,toltp(ii));
    density_ixt(ii)=(nnz(Xt)+nnz(Xt')-size(A,1))/nnz(A);
    [x5,fl5,rr5,it5,rv5] = pcg(A,b,tol,maxit,Xt,Xt',x0);
    iter_ixt(ii)=it5;
end
density_ixt
iter_ixt
semilogy(0:length(rv2)-1,rv2/norm(b),'-s',0:length(rv3)-1,rv3/norm(b),'-*',...
    0:length(rv4)-1,rv4/norm(b),'->b', 0:length(rv5)-1,rv5/norm(b),'-<r',...
    'LineWidth',1.15)
yline(tol,'r--')
plotformat(1.5,6)
legend('ichol','ixx','icholt','ixxt','Tolerance','Location','NorthEast')
xlabel('Iteration number')
ylabel('Relative residual')
if pcase==1
    title('2D Poisson equation, symmetric Legendre spectral matrix');
elseif pcase==2
    title('2D Poisson equation, finite difference matrix');
elseif pcase==3
    title('3D Poisson equation, finite difference matrix');
end








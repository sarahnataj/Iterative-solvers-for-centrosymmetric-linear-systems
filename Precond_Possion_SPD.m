function Precond_Possion_SPD(n,ecase,tol,toltp)
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
%
if ecase==0
    [A,uex]=PoissonSymLeg(n);
    elseif ecase==1
    [A,uex]=Neumann2D(n);
elseif ecase==2
    h=1/n;A=toeplitz([2; -1; zeros(n-3,1)])/(h^2);
    A=kron(speye(n-1),A)+kron(A,speye(n-1));
elseif ecase==3
    h=1/n;A=toeplitz([2; -1; zeros(n-3,1)])/(h^2);
    A=kron(kron(speye(n-1),A),speye(n-1))+...
        kron(kron(A,speye(n-1)),speye(n-1))...
        +kron(kron(speye(n-1),speye(n-1)),A);
end
x0 = zeros(size(A,1),1);
b=randn(size(A,1),1);
%b=A*uex;
maxit = size(A,1);
disp('Condition number of the matrix+++++++++++++++++++++++++++++');
condition=condest(A)
L0=chol(A,'lower');
density_chol=(nnz(L0)+nnz(L0')-size(A,1))/nnz(A)
disp('PCG without preconditioner+++++++++++++++++++++++++++++++++');
[x1,fl1,rr1,it1,rv1] = pcg(A,b,tol,maxit,[],[],x0);
it_non=it1
%
disp('PCG with ICHOL(0) factor as preconditioner++++++++++++++++++++')
L=ichol(A,struct('shape','lower'));
[x2,fl2,rr2,it2,rv2] = pcg(A,b,tol,maxit,L,L',x0);
density_ichol=(nnz(L)+nnz(L')-size(A,1))/nnz(A)
it_ichol=it2
%
disp('PCG with IXX(0) factor as preconditioner++++++++++++++++++++++')
X=ixx(A);
[x3,fl3,rr3,it3,rv3] = pcg(A,b,tol,maxit,X,X',x0);
density_ixx=(nnz(X)+nnz(X')-size(A,1))/nnz(A)
it_ixx=it3
%
disp('PCG with ICHOLT fator as preconditioner++++++++++++++++++++')
    Lt=ichol(A,struct('shape','lower','type','ict','droptol',toltp,'michol','on'));
    density_ict=(nnz(Lt)+nnz(Lt')-size(A,1))/nnz(A)
    [x4,fl4,rr4,it4,rv4] = pcg(A,b,tol,maxit,Lt,Lt',x0);
    iter_icholt=it4
%
disp('PCG with IXXT factor as preconditioner++++++++++++++++++++++')
    Xt=ixxt(A,toltp);
    density_ixt=(nnz(Xt)+nnz(Xt')-size(A,1))/nnz(A)
    [x5,fl5,rr5,it5,rv5] = pcg(A,b,tol,maxit,Xt,Xt',x0);
    iter_ixxt=it5

semilogy(0:length(rv2)-1,rv2/norm(b),'-s',0:length(rv3)-1,rv3/norm(b),'-*',...
    0:length(rv4)-1,rv4/norm(b),'->b', 0:length(rv5)-1,rv5/norm(b),'-<r',...
    'LineWidth',1.15)
yline(tol,'r--')
plotformat(1.5,6)
legend('ICHOL(0)','IXX(0)','ICHOLT(10^{-3})','IXXT(10^{-3})','tolerance','Location','NorthEast')
xlabel('Iteration number')
ylabel('Relative residual')
if ecase==0
    title('2D Poisson equation, Dirichlet boundary conditions');
elseif ecase==1
    title('2D Neumann problem');
elseif ecase==2
    title('2D Poisson equation, finite difference matrix');
elseif ecase==3
    title('3D Poisson equation, finite difference matrix');
end

end






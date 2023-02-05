function Precond_HelmholtzRobin2D(n,a, tol, toltp)
% GMRES preconditioned by incomplete double-cone factorization
% for 3D  Helmholtz equation -(\delta +a)u=f with robin boundary conditions
% Here, A is spectral differentiation matrix associated with Helmholtz equation.
% Input: n is the number of the collocation nodes
%       toltp = ilutp drop tolerance
%       tol = GMRES tolerance
%       a is square of the wave number (here we choose close to eigenvalues
%       of the 3D discrete Laplacian spectral operator with multiplicity 6
%       for the case n=20)
% Output: the number of iterations for GMRES and density of proposed preconditioners.
% Author: Sarah Nataj, email:sarah.nataj@gmail.com
%
    A=HelmholtzRobin2D(n,a);
    collocation_nodes=n+1
    size_A=size(A)
    %r=eig(full(A));min(r), max(r)
    disp('Preconditioner is centrosymmetric part of the matrix')
    condition=cond(full(A))
    J=flipud(speye(size(A,1)));
    M=0.5*(A+J*A*J);
    %s=eig(full(M\A));min(s), max(s)
    %conditionCentro(i)=cond(full(M));
    x0 = zeros(size(A,1),1);
    b=randn(size(A,1),1);
    maxit = size(A,1);
    [L, U]=lu(A);
    density_lu=(nnz(L)+nnz(U)-size(A,1))/nnz(A)
    disp('GMRES without preconditioner+++++++++++++++++++++++++++++++');
    [x0,fl0,rr0,it0,rv0] = gmres(A,b,[],tol,maxit,[],[],x0);
    iter_non=it0(2)
   disp('GMRES with ILU(0) factors as preconditioner+++++++++++++++++++');
    [L1, U1]=ilu(M);
    density_ilu=(nnz(L1)+nnz(U1)-size(A,1))/nnz(A);
    [x1,fl1,rr1,it1,rv1] = gmres(A,b,[],tol,maxit,L1,U1);
    iter_ilu=it1(2)
    %
    disp('GMRES with IXY(0) factors as preconditioner++++++++++++++++++++')
    [X1, Y1]=ixy(M);
    density_ixy=(nnz(X1)+nnz(Y1)-size(A,1))/nnz(A)
    [x2,fl2,rr2,it2,rv2]= gmres(A,b,[],tol,maxit,X1,Y1);
    iter_ixy=it2(2)
    disp('GMRES with ILUTP factors as preconditioner+++++++++++++++++');
        [L2, U2,P]=ilu(M,struct('type','ilutp','droptol',toltp));
        density_ilutp=(nnz(L2)+nnz(U2)-size(A,1))/nnz(P*A)
        [x3,fl3,rr3,it3,rv3]= gmres(P*A,P*b,[],tol,maxit,L2,U2);
        iter_ilutp=it3(2)
    disp('GMRES with IXYTP factors as preconditioner++++++++++++++++++')
        [X2, Y2, P1]=ixytp(M,toltp);
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
title('Helmholtz equation with Robin boundary conditions');

% figure(2); plot(r,0,'*r'), title('unpreconditioned')

% figure(3); plot(s,0,'*b'), title('centrosymmetric preconditioned')

plotformat(1.5,6)





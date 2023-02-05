function Precond_Biharmonic(n,k, precase, tol,toltp)
% GMRES preconditioned by incomplete double-cone factorization
% for 2D biharmonic operator and 2D biharmonic operator with variable-coefficient.
% Using equilibration algorithm to reduce the condition number.
%
% Here A is fourth order spectral differentiation (pseudospectral) matrix
% Preconditioners are:
%       precase=1; G=squared 2D second order spectral differentiation matrix (square of Laplacian),
%       precase=2; C=sparse part of the biharmonic operator,
%       precase=3; M=sparse part + dense part of A approximated by finite difference at the collocation points.
%
% Input: n is the number of the collocation nodes,
%       toltp = ilutp drop tolerance,
%       tol = GMRES tolerance,
%       k>=0 is the constant in the a(x,y)=1+kx^2y^2 for biharmonic operator with 
%       variable-coefficient, k=0 for biharmonic operator
%       precase =  preconditioner: 1 for G, 2 for C, 3 for M,
% Output: the number of iterations for GMRES and density of proposed preconditioners.
% Author: Sarah Nataj, email:sarah.nataj@gmail.com

method=2; % 1=Legendre  2=Cheb
[D, node]=PSDirv(n,method);
[A,C,M]=Biharmonic_varcoef(n,k);AA=A;
collocation_nodes=n+1
size_A=size(A)
G=Poisson_varcoef(n,2,k);
[R,S]=scaling(A);
b=randn(size(A,1),1);
b=R*b;A=R*A*S;% y=S\x the final solution.
maxit = size(A,1);
x0 = zeros(size(A,2),1);
disp('Condition number of the matrix+++++++++++++++++++++++++++++');
condition=cond(full(AA))
disp('Condition number of the matrix after equilibration++++++++++');
condition_scaled=cond(full(A))
disp('GMRES without preconditioner+++++++++++++++++++++++++++++++');
[x1,fl1,rr1,it1,rv1] = gmres(A,b,[],tol,maxit,[],[],x0);
it_non=it1(2),
%
if precase==1
    disp('preconditioner is square of Laplacian operator')
    [L01,U01]=lu(G);
    density_lu_G=(nnz(L01)+nnz(U01)-size(G,1))/nnz(G)
    disp('GMRES with ILU(0) factors of G as preconditioner+++++++++++++');
    [L,U]=ilu(G);
    density_ilu=(nnz(L)+nnz(U)-size(G,1))/nnz(G)
    [x3,fl3,rr3,it3,rv3] = gmres(A,b,[],tol,maxit,@(x)mfun(x,L,U),[],x0);
    it_ilu=it3(2),
    %
    disp('GMRES with IXY(0) factors of G as preconditioner++++++++++++++');
    [X,Y]=ixy(G);
    [x2,fl2,rr2,it2,rv2] = gmres(A,b,[],tol,maxit, @(x)mfun(x,X,Y),[],x0);
    density_ixy=(nnz(X)+nnz(Y)-size(G,1))/nnz(G)
    it_ixy=it2(2),
    %
    disp('GMRES with ILUTP factors of G as preconditioner++++++++++++');
    for ii=1:size(toltp,2)
        [L1, U1,P]=ilu(G,struct('type','ilutp','droptol',toltp(ii)));
        density_ilutp(ii)=(nnz(L1)+nnz(U1)-size(G,1))/nnz(G);
        [x3t,fl3t,rr3t,it3t,rv3t] = gmres(P*A,P*b,[],tol,maxit,@(x)mfun(x,L1,U1),[],x0);
        iter_ilutp(ii)=it3t(2);
    end
    density_ilutp
    iter_ilutp
    %
    disp('GMRES with IXYTP factors of G as preconditioner++++++++++++');
    for ii=1:size(toltp,2)
        [X1,Y1,P1]=ixytp(G,toltp(ii));
        density_ixytp(ii)=(nnz(X1)+nnz(Y1)-size(G,1))/nnz(G);
        [x2t,fl2t,rr2t,it2t,rv2t] = gmres(P1*A,P1*b,[],tol,maxit, @(x)mfun(x,X1,Y1),[],x0);
        iter_ixytp(ii)=it2t(2);
    end
    density_ixytp
    iter_ixytp
    %
    semilogy(0:length(rv3)-1,rv3/norm(U\(L\(U\(L\b)))),'-s',...
        0:length(rv2)-1,rv2/norm(Y\(X\(Y\(X\b)))),'-*',...
        0:length(rv3t)-1,rv3t/norm(U1\(L1\(U1\(L1\(P*b))))),'->b',...
        0:length(rv2t)-1,rv2t/norm(Y1\(X1\(Y1\(X1\(P1*b))))),'-<r',...
        'LineWidth',1.15);
    yline(tol,'r--');
    axis([0 120 0.5*1e-7 1])
    legend('ILU(0)','IXY(0)', 'ILUTP(10^{-3})', 'IXYTP(10^{-3})',...
        'tolerance','Location','NorthEast')
    xlabel('Iteration number')
    ylabel('Relative residual')
    title('2D biharmonic equation, preconditioner G')
    plotformat(1.5,6)
elseif precase==2
    disp('preconditioner is sparse part of the biharmonic operator')
    [L02,U02]=lu(R*C*S);
    density_lu_C=(nnz(L02)+nnz(U02)-size(C,1))/nnz(C)
    disp('GMRES with ILU(0) factors of C as preconditioner+++++++++++++');
    [LC,UC]=ilu(sparse(R*C*S));
    density_ilu=(nnz(LC)+nnz(UC)-size(C,1))/nnz(C)
    [x5,fl5,rr5,it5,rv5] = gmres(A,b,[],tol,maxit,LC,UC,x0);
    it_ilu=it5(2),
    %
    disp('GMRES with IXY(0) factors of C as preconditioner+++++++++++++');
    [XC,YC]=ixy(sparse(R*C*S));
    density_ixy=(nnz(XC)+nnz(YC)-size(C,1))/nnz(C)
    [x4,fl4,rr4,it4,rv4] = gmres(A,b,[],tol,maxit,XC,YC,x0);
    it_ixy=it4(2),
    %
    disp('GMRES with ILUTP factors of C as preconditioner+++++++++++');
    for ii=1:size(toltp,2)
        [LC1,UC1,PC]=ilu(sparse(R*C*S),struct('type','ilutp','droptol',toltp(ii)));
        density_ilutp(ii)=(nnz(LC1)+nnz(UC1)-size(C,1))/nnz(C);
        [x5t,fl5t,rr5t,it5t,rv5t] = gmres(PC*A,PC*b,[],tol,maxit,LC1,UC1,x0);
        iter_ilutp(ii)=it5t(2);
    end
    density_ilutp
    iter_ilutp
    %
    disp('GMRES with IXYTP factors of C as preconditioner+++++++++++');
    for ii=1:size(toltp,2)
        [XC1,YC1,PC1]=ixytp(sparse(R*C*S),toltp(ii));
        density_ixytp(ii)=(nnz(XC1)+nnz(YC1)-size(C,1))/nnz(C);
        [x4t,fl4t,rr4t,it4t,rv4t] = gmres(PC1*A,PC1*b,[],tol,maxit,XC1,YC1,x0);
        iter_ixytp(ii)=it4t(2);
    end
    density_ixytp
    iter_ixytp
    %
    semilogy(0:length(rv5)-1,rv5/norm(UC\(LC\b)),'-s',...
        0:length(rv4)-1,rv4/norm(YC\(XC\b)),'-*',...
        0:length(rv5t)-1,rv5t/norm(UC1\(LC1\(PC*b))),'->b',...
        0:length(rv4t)-1,rv4t/norm(YC1\(XC1\(PC1*b))),'-<r',...
        'LineWidth',1.15);
    yline(tol,'r--');
    axis([0 120 0.5*1e-7 1])
    legend('ILU(0)','IXY(0)', 'ILUTP(10^{-3})', 'IXYTP(10^{-3})', ...
        'tolerance','Location','NorthEast')
    xlabel('Iteration number')
    ylabel('Relative residual')
    title('2D biharmonic equation, preconditioner C')
    plotformat(1.5,6)
elseif precase==3
    disp('preconditioner is sparse part + dense part of A approximated by finite difference at the collocation points')
    [L03,U03]=lu(R*M*S);
    density_lu_M=(nnz(L03)+nnz(U03)-size(M,1))/nnz(M)
    disp('GMRES with ILU(0) factors of M as preconditioner++++++++++++++')
    [LM,UM]=ilu(sparse(R*M*S));
    density_ilu=(nnz(LM)+nnz(UM)-size(M,1))/nnz(M)
    [x7,fl7,rr7,it7,rv7] = gmres(A,b,[],tol,maxit,LM,UM,x0);
    it_ilu=it7(2),
    %
    disp('GMRES with IXY(0) factors of M as preconditioner++++++++++++++')
    [XM,YM]=ixy(sparse(R*M*S));
    density_ixy=(nnz(XM)+nnz(YM)-size(M,1))/nnz(M)
    [x6,fl6,rr6,it6,rv6] = gmres(A,b,[],tol,maxit,XM,YM,x0);
    it_ixy=it6(2),
    %
    disp('GMRES with ILUTP factors of M as preconditioner++++++++++++')
    for ii=1:size(toltp,2)
        [LM1,UM1,PM]=ilu(sparse(R*M*S),struct('type','ilutp','droptol',toltp(ii)));
        density_ilutp(ii)=(nnz(LM1)+nnz(UM1)-size(M,1))/nnz(M);
        [x7t,fl7t,rr7t,it7t,rv7t] = gmres(PM*A,PM*b,[],tol,maxit,LM1,UM1,x0);
        iter_ilutp(ii)=it7t(2);
    end
    density_ilutp
    iter_ilutp
    %
    disp('GMRES with IXYTP factors of M as preconditioner++++++++++++')
    for ii=1:size(toltp,2)
        [XM1,YM1,PM1]=ixytp(sparse(R*M*S),toltp(ii));
        density_ixytp(ii)=(nnz(XM1)+nnz(YM1)-size(M,1))/nnz(M);
        [x6t,fl6t,rr6t,it6t,rv6t] = gmres(PM1*A,PM1*b,[],tol,maxit,XM1,YM1,x0);
        iter_ixytp(ii)=it6t(2);
    end
    density_ixytp
    iter_ixytp
    %
    semilogy(0:length(rv7)-1,rv7/norm(UM\(LM\b)),'-s',...
        0:length(rv6)-1,rv6/norm(YM\(XM\b)),'-*',...
        0:length(rv7t)-1,rv7t/norm(UM1\(LM1\(PM*b))),'->b',...
        0:length(rv6t)-1,rv6t/norm(YM1\(XM1\(PM1*b))),'-<r',...
        'LineWidth',1.15);
    yline(tol,'r--');
    axis([0 120 0.5*1e-7 1])
    legend('ILU(0)','IXY(0)', 'ILUTP(10^{-3})','IXYTP(10^{-3})',...
        'tolerance','Location','NorthEast')
    xlabel('Iteration number')
    ylabel('Relative residual')
    title('2D biharmonic equation, preconditioner M')
    plotformat(1.5,6)
end
end
function M=fdpoisson(node)
n=size(node,1)-1;
for j=1:size(node)-1
    h(j)=node(j)-node(j+1);
end
B=zeros(n-1,n-1);
B(1,1)=2/(h(2)*h(1));
B(1,2)=-2/(h(2)*(h(2)+h(1)));
for j=2:n-2
    B(j,j)=2/(h(j)*h(j+1));
    B(j,j+1)=-2/(h(j+1)*(h(j)+h(j+1)));
    B(j,j-1)=-2/(h(j)*(h(j)+h(j+1)));
end
B(n-1,n-2)=-2/(h(n-1)*(h(n-1)+h(n)));
B(n-1,n-1)=2/(h(n)*h(n-1));
M=B;
end

function y=mfun(b,L,U)
y=U\(L\(U\(L\b)));
end





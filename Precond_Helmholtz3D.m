function Precond_Helmholtz3D(n,a,tol)
% GMRES preconditioned by incomplete double-cone factorization
% for 3D  Helmholtz equation -(\Delta +a)u=f
% Here, A is  spectral differentiation matrix associated with Helmholtz equation.
% Input: n is the number of the collocation nodes
%       toltp = ilutp drop tolerance
%       tol = GMRES tolerance
%       a is square of the wave number
% Output: the number of iterations for GMRES and density of proposed preconditioners.
% Author: Sarah Nataj, email:sarah.nataj@gmail.com

[A, b]=Helmholtz3D(n+1,a);
%condition=cond(full(A))
%size(full(A),1)
wave_number=a
collecation_nodes=n+2
size_A=size(A)
condition_A=cond(full(A))
x00 = zeros(size(A,1),1);
%x00 = randn(size(A,1),1);
%b=randn(size(A,1),1);
maxit=200;
re=25;
%maxit_re = ceil(size(A,1)/re);
maxit_re = ceil(maxit/re);
%[L, U]=lu(A);
%density_lu=(nnz(L)+nnz(U)-size(A,1))/nnz(A)

%disp('GMRES without preconditioner+++++++++++++++++++++++++++++++');
l=1;
%tic
%[x0,fl0,rr0,it0,rv0] = gmres(A,b,[],tol,maxit);
%time(l)=toc;
%iter(l)=it0(2);

%disp('ReStart');
%tic
%[x0,fl0,rr0,it0,rv0] = gmres(A,b,re,tol,maxit_re);
%time_re(l)=toc;
%iter_re(l)=(it0(1)-1)*re+it0(2);
%l=l+1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%disp('GMRES with ILU(0) factors as preconditioner+++++++++++++++++++');
[L1, U1]=ilu(A);
%density_ilu=(nnz(L1)+nnz(U1)-size(A,1))/nnz(A)
tic
[x1,fl1,rr1,it1,rv1] = gmres(A,b,[],tol,maxit,L1,U1,x00);
time(l)=toc;
iter(l)=it1(2);

tic
[x1a,fl1a,rr1a,it1a,rv1a] = gmres(A,b,re,tol,maxit_re,L1,U1,x00);
time_re(l)=toc;
iter_re(l)=(it1a(1)-1)*re+it1a(2);
l=l+1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%disp('GMRES with IXY(0) factors as preconditioner++++++++++++++++++++')
[X1, Y1]=ixy(A);
%density_ixy=(nnz(X1)+nnz(Y1)-size(A,1))/nnz(A)
tic
[x2,fl2,rr2,it2,rv2]= gmres(A,b,[],tol,maxit,X1,Y1,x00);
time(l)=toc;
iter(l)=it2(2);

%disp('ReStart')
tic
[x2a,fl2a,rr2a,it2a,rv2a]= gmres(A,b,re,tol,maxit_re,X1,Y1,x00);
time_re(l)=toc;
iter_re(l)=(it2a(1)-1)*re+it2a(2);
l=l+1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% toltp=1e-3;
% %disp('GMRES with ILUTP factors as preconditioner+++++++++++++++++');
% [L2, U2,P]=ilu(A,struct('type','ilutp','droptol',toltp));
% %density_ilutp=(nnz(L2)+nnz(U2)-size(A,1))/nnz(P*A)
% tic
% [x3,fl3,rr3,it3,rv3]= gmres(P*A,P*b,[],tol,maxit,L2,U2,x00);
% time(l)=toc;
% iter(l)=it3(2);
%
% %disp('ReStart');
% [L2, U2,P]=ilu(A,struct('type','ilutp','droptol',toltp));
% tic
% [x3,fl3,rr3,it3,rv3]= gmres(P*A,P*b,re,tol,maxit_re,L2,U2,x00);
% time_re(l)=toc;
% iter_re(l)=(it3(1)-1)*re+it3(2);
% l=l+1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%disp('GMRES with IXYTP factors as preconditioner++++++++++++++++++')
%[X2, Y2, P1]=ixytp(A,toltp);
%density_ixytp=(nnz(X2)+nnz(Y2)-size(A,1))/nnz(P1*A)
%tic
%[x4,fl4,rr4,it4,rv4]= gmres(P1*A,P1*b,[],tol,maxit,X2,Y2,x00);
%time(l)=toc;
%iter(l)=it4(2);

%disp('ReStart')
%tic
%[x4a,fl4a,rr4a,it4a,rv4a]= gmres(P1*A,P1*b,re,tol,maxit_re,X2,Y2,x00);
%time_re(l)=toc;
%iter_re(l)=(it4(1)-1)*re+it4(2);
%l=l+1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
toltp=1e-4;
%disp('GMRES with ILUTP factors as preconditioner+++++++++++++++++');
[L2, U2,P]=ilu(A,struct('type','ilutp','droptol',toltp));
%density_ilutp=(nnz(L2)+nnz(U2)-size(A,1))/nnz(P*A)
tic
[x3,fl3,rr3,it3,rv3]= gmres(P*A,P*b,[],tol,maxit,L2,U2,x00);
time(l)=toc;
iter(l)=it3(2);

%disp('ReStart');
[L2, U2,P]=ilu(A,struct('type','ilutp','droptol',toltp));
tic
[x3a,fl3a,rr3a,it3a,rv3a]= gmres(P*A,P*b,re,tol,maxit_re,L2,U2,x00);
time_re(l)=toc;
iter_re(l)=(it3a(1)-1)*re+it3a(2);
l=l+1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%disp('GMRES with IXYTP factors as preconditioner++++++++++++++++++')
[X2, Y2, P1]=ixytp(A,toltp);
%density_ixytp=(nnz(X2)+nnz(Y2)-size(A,1))/nnz(P1*A)
tic
[x4,fl4,rr4,it4,rv4]= gmres(P1*A,P1*b,[],tol,maxit,X2,Y2,x00);
time(l)=toc;
iter(l)=it4(2);

%disp('ReStart')
tic
[x4a,fl4a,rr4a,it4a,rv4a]= gmres(P1*A,P1*b,re,tol,maxit_re,X2,Y2,x00);
time_re(l)=toc;
iter_re(l)=(it4a(1)-1)*re+it4a(2);


figure (1); semilogy(0:length(rv1)-1,rv1/norm(U1\(L1\b)),'-s',...
    0:length(rv2)-1,rv2/norm(Y1\(X1\b)),'-*',...
    0:length(rv3)-1,rv3/norm(U2\(L2\b)),'->b',...
    0:length(rv4)-1,rv4/norm(Y2\(X2\b)),'-<r',...
    'LineWidth',1.15)
yline(tol,'r--')
plotformat(1.5,6)
legend('ILU(0)','IXY(0)','ILUTP(10^{-4})','IXYTP(10^{-4})','tolerance','Location','NorthEast')
xlabel('Iteration number')
ylabel('Relative residual')
title('3D Helmholtz equation with Dirichlet boundary conditions');
plotformat(1.5,6)

figure (2); semilogy(0:length(rv1a)-1,rv1a/norm(U1\(L1\b)),'-s',...
    0:length(rv2a)-1,rv2a/norm(Y1\(X1\b)),'-*',...
    0:length(rv3a)-1,rv3a/norm(U2\(L2\b)),'->b',...
    0:length(rv4a)-1,rv4a/norm(Y2\(X2\b)),'-<r',...
    'LineWidth',1.15)
yline(tol,'r--')
plotformat(1.5,6)
legend('ILU(0)','IXY(0)','ILUTP(10^{-4})','IXYTP(10^{-4})','tolerance','Location','NorthEast')
xlabel('Iteration number')
ylabel('Relative residual')
title('Restarted 3D Helmholtz equation with Dirichlet boundary conditions');
plotformat(1.5,6)

Table=table(iter',iter_re',time',time_re','VariableNames',{'Iter no restart', 'Iter with restart','Time no restart', 'Time with restart'},...
    RowNames={'ILU(0)', 'IXY(0)','ILUTP(10^-4)','IXYTP(10^-4)'})
end


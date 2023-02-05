function [Nnop,Nilu,Nixy,Nilutp,Nixytp]=gmres_biharmonic(n,k, precase)
% GMRES preconditioned by incomplete double-cone factorization
% for 2D biharmonic operator and 2D biharmonic operator with variable-coefficient.
% Use equilibration algorithm to reduce the condition number.
%
% Here A is fourth order spectral differentiation (pseudospectral) matrix
% Preconditioners are:
%       precase=1; G=squared 2D second order spectral differentiation matrix (square of Laplacian),
%       precase=2; C=sparse part of the biharmonic operator,
%       precase=3; M=sparse part + dense part of A approximated by finite difference at the collocation points.
%
%Input: n is the number of the collocation nodes,
%       toltp = ilutp drop tolerance,
%       tol = GMRES tolerance,
%       k>0 is the constant in the a(x,y) for biharmonic operator with
%       variable-coefficient, k=0 for biharmonic operator
%       precase =  preconditioner: 1 for G, 2 for C, 3 for M,
%Output: the number of iterations for GMRES and density of proposed preconditioners.
% Author: Sarah Nataj, email:sarah.nataj@gmail.com
tol = 1e-6;
toltp=1e-3;
[A,C,M]=Biharmonic_varcoef(n,k);
[R,S]=scaling(A);
b=randn(size(A,1),1);
b=R*b;A=R*A*S;% y=S\x the final solution.
maxit = size(A,1);
x0 = zeros(size(A,2),1);
[x1,fl1,rr1,it1,rv1] = gmres(A,b,[],tol,maxit,[],[],x0);
%size(full(A),1)
Nnop=it1(2);
%
if precase==2
    [LC,UC]=ilu(sparse(R*C*S));
    [x5,fl5,rr5,it5,rv5] = gmres(A,b,[],tol,maxit,LC,UC,x0);
    Nilu=it5(2);
    %
    [XC,YC]=ixy(sparse(R*C*S));
    [x4,fl4,rr4,it4,rv4] = gmres(A,b,[],tol,maxit,XC,YC,x0);
    Nixy=it4(2);
    %
    [LC1,UC1,PC]=ilu(sparse(R*C*S),struct('type','ilutp','droptol',toltp));
    [x5t,fl5t,rr5t,it5t,rv5t] = gmres(PC*A,PC*b,[],tol,maxit,LC1,UC1,x0);
    Nilutp=it5t(2);
    %
    [XC1,YC1,PC1]=ixytp(sparse(R*C*S),toltp);
    [x4t,fl4t,rr4t,it4t,rv4t] = gmres(PC1*A,PC1*b,[],tol,maxit,XC1,YC1,x0);
    Nixytp=it4t(2);

elseif precase==3
    [LM,UM]=ilu(sparse(R*M*S));
    [x7,fl7,rr7,it7,rv7] = gmres(A,b,[],tol,maxit,LM,UM,x0);
    Nilu=it7(2);
    %
    [XM,YM]=ixy(sparse(R*M*S));
    [x6,fl6,rr6,it6,rv6] = gmres(A,b,[],tol,maxit,XM,YM,x0);
    Nixy=it6(2);
    %
    [LM1,UM1,PM]=ilu(sparse(R*M*S),struct('type','ilutp','droptol',toltp));
    [x7t,fl7t,rr7t,it7t,rv7t] = gmres(PM*A,PM*b,[],tol,maxit,LM1,UM1,x0);
    Nilutp=it7t(2);
    %
    [XM1,YM1,PM1]=ixytp(sparse(R*M*S),toltp);
    [x6t,fl6t,rr6t,it6t,rv6t] = gmres(PM1*A,PM1*b,[],tol,maxit,XM1,YM1,x0);
    Nixytp=it6t(2);
end






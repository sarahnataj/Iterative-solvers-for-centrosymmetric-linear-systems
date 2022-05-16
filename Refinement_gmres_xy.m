% Given A, b in be given in precision u, this program solves 
% Ax = b using GMRES-based iterative refinement with inexact XY factorization 
% (calculated in single precision).( Algorithm 5).
%
% Input: A is the spectral differentiation (pseudospectral) matrix of size n
%       associated with
%            1D Poisson/variable--coefficient Poisson (diffusion) equations,
%            1D/2D biharmonic/variable--coefficient biharmonic equations,
%            3D Helmholtz equation.
%       u = the working precision at which A and b are stored=double,
%       uf = the precision at which the factorization of 𝐴 is computed = single,
%       ur = the precision at which the residual is calculated = double,
%       tol is the convergence criterion for the refinement process,
%       tolg is the tolerance for convergence of GMRES.
% Output: n_ir = the number of iterations for GMRES to converge,
%       e_u = norm of the error.
%
% Algorithm 5 is discussed in the manuscript: Structure-preserving 
% solvers for centrosymmetric linear systems with applications to spectral methods
% by Chen Greif, Sarah Nataj and Manfred Trummer.
%
% Author: Sarah Nataj, email:sarah.nataj@gmail.com
clc
clear all;
close all; 
tol=1e-8;
tolg=1e-3;
scale=1;% use scaling to reduce the condition number
imax=100;
%n=15;[A, uex]=PseudoSpectral1D(n,2);%1D poisson
%n=21;[A,uex]=PseudoSpectral1D(n,4);%1D bihaarmonic
n=20;[A,uex]=Poisson_varcoef(n,1,10);% 1D variable-coefficient Poisson
%n=19;[A,uex]=PseudoSpectral2D(n,4);%2D biharmonic
%n=23;[A, C, M, uex]=Biharmonic_varcoef(n,1);A=full(A);%2D variable coefficient biharmonic
%n=15;a=600;[A,uex]=Helmholtz3D(n,a);A=full(A);
condition=cond(A)
x_ext=uex;
b=A*x_ext;
maxit = size(A,1);
if scale==1
    [R,S]=scaling(A);
    Ah=R*A*S;bh=R*b;
    condition_scaled=cond(Ah)
else
    Ah=A;bh=b;
end
[P,X,Y]=xy(single(Ah),0);
w=SolveConeX(X,P*bh);
y_uf=SolveConeY(Y,w);
if scale==1
x_u=double(S*y_uf);
else
 x_u=double(y_uf);
end
e_ur = double(x_ext-x_u);
disp('Error from solving by inexact XY factorization+++++++++')
e_u = double(norm(e_ur))
disp('Error after using iterative refinement+++++++++++++++++')
% solve the residual equation
for i=1:imax
r_u = double(b-A*x_u);
if scale==1
    r_u=R*r_u;
    [z_u,fl,rr,it,rv]=gmres(P*Ah,P*r_u,[],tolg,maxit,X,Y);z_u=S*z_u;
else
    [z_u,fl,rr,it,rv]=gmres(P*A,P*r_u,[],tolg,maxit,X,Y);
end
x_u=x_u+z_u;
e_ur = double(x_ext-x_u);
e_u = double(norm(e_ur));
if e_u<=tol
    x_u;e_u,n_ir=it(2)
    return 
end

end

function x=SolveConeX(X,b)
X=single(X);b=single(b);
n=size(X,1);
x=single(zeros(n,1));
m=floor(n/2);
for k=1:m
    p=k;
    q=n-k+1;
    XX=[X(p,p) X(p,q);X(p,q) X(p,p)];
    if (p==1)
        bb=[b(p)  b(q)]';
    else
        bb=[b(p)-X(p,1:p-1)*x(1:p-1)-X(p,q+1:end)*x(q+1:end),...
            b(q)-X(q,1:p-1)*x(1:p-1)-X(q,q+1:end)*x(q+1:end)]';
    end
    xx=XX\bb;
    x(p)=xx(1);x(q)=xx(2);
end
if mod(n,2)==1
    x(m+1)=single((b(m+1)-X(m+1,:)*x)/X(m+1,m+1));
end
end

function x=SolveConeY(X,b)
X=single(X);b=single(b);
n=size(X,1);
x=single(zeros(n,1));
m=floor(n/2);
if mod(n,2)==1
    x(m+1)=b(m+1)/X(m+1,m+1);
end
for k=1:m
    p=m-k+1;
    q=m+k+mod(n,2);
    XX=[X(p,p) X(p,q);X(p,q) X(p,p)];
    if (p==m)&&(mod(n,2)==0)
        bb=[b(p)  b(q)]';
    else
        bb=[b(p)-X(p,p+1:q-1)*x(p+1:q-1)  b(q)-X(q,p+1:q-1)*x(p+1:q-1)]';
    end
    xx=XX\bb;
    x(p)=xx(1);x(q)=xx(2);
end
end






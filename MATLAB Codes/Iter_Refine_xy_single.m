% Given A, b, and x_ext in be given in precision u. This algorithm solves 
% Ax = b using GMRES-based iterative refinement using XY factorization 
% calculated in single precision.( Algorithm 5).
%
% Input: A is the spectral differentiation (pseudospectral) matrix of size n
%       associated with
%            1D Poisson/variable--coefficient Poisson (diffusion),
%            1D biharmonic/variable--coefficient biharmonic,
%            2D biharmonic.
%       b is given right hand side.
%       x_ext= the exact solution (here we assume a random vector of size n)
%       u = the working precision at which A and b are stored=double,
%       uf = the precision at which the factorization of 𝐴 is computed = single,
%       ur = the precision at which the residual is calculated = double,
%       tol is the convergence criterion for the refinement process,
%       tolg is the tolerance for convergence of GMRES.
% Output: it = the number of iterations for GMRES to converge,
%       e_u = norm of the error.
% Author: Sarah Nataj, email:sarah.nataj@gmail.com
clc
clear all;
close all; 
tol=1e-8;
tolg=1e-4;
imax=100;
n=20;
%A=PseudoSpectral1D(n+1,2);%1D poisson
%A=PseudoSpectral1D(n+1,4);%1D bihaarmonic
%A=PseudoSpectral_Poisson_varcoef(n+1,1,10);% 1D variable-coefficient Poisson
%A=PseudoSpectral2D(n+1,4);%2D biharmonic
[A, C, M]=PseudoSpectral_Biharmonic_varcoef(n+1,10);A=full(A);%2D variable coefficient biharmonic
x_ext=rand(size(A,1),1);
b=A*x_ext;
maxit = size(A,1);
[R,S]=scaling(A);
Ah=R*A*S;bh=R*b;
[P,X,Y]=xy(single(Ah),0);
w=SolveConeX(X,P*bh);
y_uf=SolveConeY(Y,w);
x_u=double(S*y_uf);
e_ur = double(x_ext-x_u);
e_u = double(norm(e_ur));
% solve the residual equation
B=S*(Y\(X\(P*R*A)));BB=S*(Y\(X\(P*R)));
for i=1:imax
r_u = double(b-A*x_u);
bb=BB*r_u;
[z_u,fl,rr,it,rv]=gmres(B,bb,[],tolg,maxit);
x_u=x_u+z_u;
e_ur = double(x_ext-x_u);
e_u = double(norm(e_ur));
if e_u<=tol
    x_u;e_u,it,
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






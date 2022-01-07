% This program solves a nonsingular square linear sytem Az=b of size n by n
% for a symmetric positive definite (SPD) centrosymmetric matrix A using 
% double-cone factorization and modified backward substitution (Algorithm 2).
%
% Here A is 2D second order SPD spectral differentiation matrix using
% Legendre Lobatto points. We also assume SPD matrices arising from solving
% 2D/3D Poisson equation using finite differences schemes. The exact solution 
% z_exact is the vector with all ones. 
%
% Input: nn is the array of number of collocation nodes to be considered.
% Output: Graphs for time complexity of the algorithm and relative errors 
%        with respect to n, the size of the matrix.
%
% For large n, change xx(A,1) to xx(A,0) in line 38 and pchol to chol in line 45.
% Author: Sarah Nataj, email:sarah.nataj@gmail.com
clc
clear all
close all
nn=[10:10:30];ecase=1;%symmetric Legendre spectral matrix for 2D Poisson
%nn=[10:10:60];ecase=2;%finite difference for 2D Poisson
%nn=[3:2:9];ecase=3;%finite difference for 3D Poisson
for ii=1:size(nn,2)
    if ecase==1
    A=SymLeg(nn(ii)+1);A=full(A);
    elseif ecase==2
    n=nn(ii);h=1/n;A=toeplitz([2; -1; zeros(n-3,1)])/(h^2);
    A=kron(eye(n-1),A)+kron(A,eye(n-1));
    elseif ecase==3
    n=nn(ii);h=1/n;A=toeplitz([2; -1; zeros(n-3,1)])/(h^2);
    A=kron(kron(eye(n-1),A),eye(n-1))+kron(kron(A,eye(n-1)),eye(n-1))+kron(kron(eye(n-1),eye(n-1)),A);
    end
    n=size(A,1);A=0.5*(A+A');J=flipud(eye(n));
    A=0.5*(A+J*A*J);
    n=size(A,1);nn(ii)=n;
    z_exact=ones(n,1);
    b=A*z_exact;AA=A;bb=b;
    tic
    X=xx(A,1);%X=xx(A,0);
    %Solving Az=b, since XX'z=b, therfore solve Xw=b, and then X'z=w
    w=SolveConeX(X,b);
    z=SolveConeY(X',w);
    time1(ii)=toc;
    norm1(ii)=norm(z_exact-z)/norm(z_exact);
    tic
    RR=pchol(AA);%RR=chol(AA);
    ww=ForsubL(RR',bb);
    zz=BacksubU(RR,ww);
    time2(ii)=toc;
    norm2(ii)=norm(z_exact-zz)/norm(z_exact);  
end

figure(1); plot(nn,time1,'-o',nn,time2,'-x','LineWidth',2)
legend('Algorithm 2','Cholesky solver','Location','NorthWest'),
titlefun(ecase);
xlabel('n')
ylabel('Time')
plotformat(1.5,6)

figure(2); loglog(nn,norm1,'-o',nn,norm2,'-x','LineWidth',2)
legend('Algorithm 2','Cholesky solver','Location','NorthWest'),
titlefun(ecase);
ylabel('Relative error');
xlabel('n')
plotformat(1.5,6)


function L=pchol(A)
n=size(A,1);
for k=1:n-1
    A(k,k)=sqrt(A(k,k));
    for i=k+1:n
        A(i,k)=A(i,k)/A(k,k);
    end
    for j=k+1:n
        for i=j:n
            A(i,j)=A(i,j)-A(i,k)*A(j,k);
        end
    end
end
A(n,n)=sqrt(A(n,n));
L=zeros(n,n);
for k=1:n
    L(k,1:k)=A(k,1:k);
end
L=L';
end

function x=SolveConeX(X,b)
n=size(X,1);
x=zeros(n,1);
m=floor(n/2);
for k=1:m
    p=k;
    q=n-k+1;
    XX=[X(p,p) X(p,q);X(p,q) X(p,p)];
    if (p==1)
        bb=[b(p)  b(q)]';
    else
        bb=[b(p)-X(p,1:p-1)*x(1:p-1)-X(p,q+1:end)*x(q+1:end)  b(q)-X(q,1:p-1)*x(1:p-1)-X(q,q+1:end)*x(q+1:end)]';
    end
    xx=XX\bb;
    x(p)=xx(1);x(q)=xx(2);
end
if mod(n,2)==1
    x(m+1)=(b(m+1)-X(m+1,:)*x)/X(m+1,m+1);
end
end

function x=SolveConeY(X,b)
n=size(X,1);
x=zeros(n,1);
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

function x=ForsubL(A,b)
n=size(A,1);
x=zeros(n,1);
x(1)=b(1)/A(1,1);
for k=2:n
    x(k)=(b(k)-A(k,1:k-1)*x(1:k-1))/A(k,k);
end
end

function x=BacksubU(A,b)
n=size(A,1);
x=zeros(n,1);
x(n)=b(n)/A(n,n);
for k=n-1:-1:1
    x(k)=(b(k)-A(k,k+1:n)*x(k+1:n))/A(k,k);
end
end

function titlefun(ecase)
if ecase==1
    title('2D second-order symmetric Legendre spectral matrices');
elseif ecase==2
    title('Finite difference matrices for 2D Poisson equation');
elseif ecase==3
    title('Finite difference matrices for 3D Poisson equation');
end
end




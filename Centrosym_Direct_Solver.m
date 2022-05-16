% This program solves a nonsingular square linear sytem Az=b of size n by n
% for a centrosymmetric matrix A using double-cone factorization and
% modified backward substitution (Algorithm 1 and Algorithm 2).
%
% Algorithm 1: A is 2D/3D second order and 2D fourth order spectral
% differentiation matrices. z_exact=uex is exact solution of 2D/3D Poisson
% and  2D biharmonic equations at the collocation points.
%
% Algorithm 2: A is 2D second order SPD spectral differentiation matrix using
% Legendre Gauss Lobatto points. We also assume SPD matrices arising from solving
% 3D Poisson equation using finite differences schemes. 
%
% Input: N is the array of number of collocation nodes to be considered.
% Output: Graphs for time complexity of the algorithm and relative errors
%        with respect to n, the size of the matrix
%
% The Algorithms 1 and 2 are discussed in the manuscript: Structure-preserving 
% solvers for centrosymmetric linear systems with applications to spectral methods
% by Chen Greif, Sarah Nataj and Manfred Trummer.
%
% Author: Sarah Nataj, email:sarah.nataj@gmail.com
clc
clear all
close all

N=10:1:20;ecase=1;%2Dpoisson
N=10:1:20;ecase=2;%2Dbiharmonic
%N=5:1:10;ecase=3;%3Dpoisson

%N=8:2:21;ecase=10;%symmetric Legendre spectral matrix for 2D Poisson equaion with homogeneous Dirichlet boundary conditions
%N=8:2:21;ecase=11; %symmetric Legendre spectral matrix for the Neumann  problem -\laplacian u+u=f(x,y)
N=30:2:40;ecase=12;%finite difference approach for 2D Poisson equation with homogeneous Dirichlet boundary conditions
%N=5:2:13;ecase=13;%finite difference approach for 3D Poisson equation with homogeneous Dirichlet boundary conditions

if ecase<=9
    for i=1:size(N,2)
        n=N(i);
        if ecase==1
            [A,uex]=PseudoSpectral2D(n+1,2);
        elseif ecase==2
            [A,uex]=PseudoSpectral2D(n+1,4);
        elseif ecase==3
            [A,uex]=PseudoSpectral3D(n+1);
        end
        n=size(full(A),1);N(i)=n;
        z_exact=uex;
        b=A*z_exact;
        % using equilibration algorithm to improve the condition number
        [R,S]=scaling(A);A=R*A*S;b=R*b;
        % calculate the time for solving the linear system using double-cone factorization
        tic
        [P,X,Y]=xy(A,1);%[P,X,Y]=xy(A,0);
        % solving the linear system Az=b, since P'XYz=b , therefore XYz=Pb , then solve Xw=b, and then Yz=w.
        w=SolveConeX(X,P*b);
        z=SolveConeY(Y,w);
        time1(i)=toc;
        z=S*z;
        norm1(i)=norm(z_exact-z)/norm(z_exact);
        % calculate the time for solving the linear system using LU solver
        tic
        [LL,UU,PP]=plu(A);%[LL,UU,PP]=lu(A);
        ww=ForsubL(LL,PP*b);
        zz=BacksubU(UU,ww);
        time2(i)=toc;
        zz=S*zz;
        norm2(i)=norm(z_exact-zz)/norm(z_exact);
    end
    figure(1); loglog(N,norm1,'-o',N,norm2,'-x','LineWidth',2)
    legend('Algorithm 1','LU solver','Location','NorthWest'),
    titlefun1(ecase);
    xlabel('n')
    ylabel('Relative error')
    plotformat(1.5,6)

    figure(2); plot(N,time1,'-o',N,time2,'-x','LineWidth',2)
    legend('Algorithm 1','LU solver','Location','NorthWest'),
    titlefun1(ecase);
    xlabel('n')
    ylabel('Time')
    plotformat(1.5,6)

elseif ecase>=10
    for i=1:size(N,2)
        if ecase==10
            [A, uex]=PoissonSymLeg(N(i)+1);%A=full(A);
        elseif ecase==11
            [A, uex]=Neumann2D(N(i)+1);%A=full(A);
        elseif ecase==12
            n=N(i);h=1/n;A=toeplitz([2; -1; zeros(n-3,1)])/(h^2);
            A=kron(eye(n-1),A)+kron(A,eye(n-1));
            uex=ones(size(A,1),1);
        elseif ecase==13
            n=N(i);h=1/n;A=toeplitz([2; -1; zeros(n-3,1)])/(h^2);
            A=kron(kron(eye(n-1),A),eye(n-1))+kron(kron(A,eye(n-1)),eye(n-1))+kron(kron(eye(n-1),eye(n-1)),A);
            uex=ones(size(A,1),1);
        end
        n=size(A,1);N(i)=n;
        b=A*uex;
        tic
        X=xx(A,1);%X=xx(A,0);
        %Solving Az=b, since XX'z=b, therfore solve Xw=b, and then X'z=w
        w=SolveConeX(X,b);
        z=SolveConeY(X',w);
        time1(i)=toc;
        norm1(i)=norm(uex-z)/norm(uex);
        tic
        RR=pchol(A);%RR=chol(A);
        ww=ForsubL(RR',b);
        zz=BacksubU(RR,ww);
        time2(i)=toc;
        norm2(i)=norm(uex-zz)/norm(uex);
    end

    figure(1); loglog(N,norm1,'-o',N,norm2,'-x','LineWidth',2)
    legend('Algorithm 2','Cholesky solver','Location','NorthWest'),
    titlefun2(ecase);
    ylabel('Relative error');
    xlabel('n')
    plotformat(1.5,6)

    figure(2); plot(N,time1,'-o',N,time2,'-x','LineWidth',2)
    legend('Algorithm 2','Cholesky solver','Location','NorthWest'),
    titlefun2(ecase);
    xlabel('n')
    ylabel('Time')
    plotformat(1.5,6)

end

function [L, U, P]=plu(A)
n=size(A,1);
p=1:n;
for k=1:n-1
    [val,q]=max(abs(A(k:n,k)));
    q=q+(k-1);
    A([k,q],:)=A([q,k],:);
    p([k,q])=p([q,k]);
    J=k+1:n;
    A(J,k)=A(J,k)/A(k,k);
    A(J,J)=A(J,J)-A(J,k)*A(k,J);
end
L=zeros(n,n);U=zeros(n,n);
for k=1:n
    L(k,1:k-1)=A(k,1:k-1);L(k,k)=1;
    U(k,k:n)=A(k,k:n);
    I=eye(n);
    P=I(p,:);
end
end

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

function titlefun1(ecase)
if ecase==1
    title('2D Poisson equation');
elseif ecase==2
    title('2D biharmonic equation');
elseif ecase==3
    title('3D Poisson equation');
end
end

function titlefun2(ecase)
if ecase==10
    title('Dirichlet boundary conditions');
elseif    ecase==11
    title('Neumann boundary conditions');
elseif ecase==12
    title('Finite difference matrices for 2D Poisson equation');
elseif ecase==13
    title('Finite difference matrices for 3D Poisson equation');
end
end

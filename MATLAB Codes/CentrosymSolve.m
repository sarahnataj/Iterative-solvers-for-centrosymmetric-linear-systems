% This program solves a nonsingular square linear sytem Az=b of size n by n
% for a centrosymmetric matrix A using double-cone factorization and
% modified backward substitution (Algorithm 1).
%
% Here, A is 1D/2D/3D second order and 1D/2D fourth order spectral 
% differentiation matrices. z_exact=uex is exact solution of 1D/2D/3D Poisson 
% and  1D/2D biharmonic equations at the collocation points.
%
% Input: nn is the array of number of collocation nodes to be considered.
% Output: Graphs for time complexity of the algorithm and relative errors 
%        with respect to n, the size of the matrix
%
% For large n, change xy(A,1) to xy(A,0) in line 43 and plu to lu in line 53.
% Author: Sarah Nataj, email:sarah.nataj@gmail.com
clc
clear all
close all
%nn=[100:10:200];ecase=1;%1Dpoisson
%nn=[60:20:400];ecase=2;%1Dbiharmonic
%nn=[10:5:20];ecase=3;%2Dpoisson
%nn=[10:5:20];ecase=4;%2Dbiharmonic
nn=[4:2:10];ecase=5;%3Dpoisson

for ii=1:size(nn,2)
    n=nn(ii);
    if ecase==1
        [A,uex]=PseudoSpectral1D(n+1,2);
    elseif ecase==2
        [A,uex]=PseudoSpectral1D(n+1,4);
    elseif ecase==3
        [A,uex]=PseudoSpectral2D(n+1,2);%A=full(A);
    elseif ecase==4
        [A,uex]=PseudoSpectral2D(n+1,4);
    elseif ecase==5
        [A,uex]=PseudoSpectral3D(n+1);%A=full(A);
    end
    n=size(full(A),1);nn(ii)=n;
    z_exact=uex;
    b=A*z_exact;
    % using equilibration algorithm to improve the condition number
    [R,S]=scaling(A);A=R*A*S;b=R*b;bb=b;
    % calculate the time for solving the linear system using double-cone factorization
    tic
    [P,X,Y]=xy(A,1);%[P,X,Y]=xy(A,0);
    % solving the linear system Az=b, since P'XYz=b , therefore XYz=Pb , then solve Xw=b, and then Yz=w.
    b=P*b;
    w=SolveConeX(X,b);
    z=SolveConeY(Y,w);
    time1(ii)=toc;
    z=S*z;
    norm1(ii)=norm(z_exact-z)/norm(z_exact);
    % calculate the time for solving the linear system using LU solver
    tic
    [LL,UU,PP]=plu(A);%[LL,UU,PP]=lu(AA);
    ww=ForsubL(LL,PP*bb);
    zz=BacksubU(UU,ww);
    time2(ii)=toc;
    zz=S*zz;
    norm2(ii)=norm(z_exact-zz)/norm(z_exact);
end


figure(1); plot(nn,time1,'-o',nn,time2,'-x','LineWidth',2)
legend('Algorithm 1','LU solver','Location','NorthWest'),
titlefun(ecase);
xlabel('n')
ylabel('Time')
plotformat(1.5,6)



figure(2); loglog(nn,norm1,'-o',nn,norm2,'-x','LineWidth',2)
legend('Algorithm 1','LU solver','Location','NorthWest'),
titlefun(ecase);
xlabel('n')
ylabel('Relative error')
plotformat(1.5,6)


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
    title('1D Poisson equation');
elseif ecase==2
    title('1D biharmonic equation');
elseif ecase==3
    title('2D Poisson equation');
elseif ecase==4
    title('2D biharmonic equation');
elseif ecase==5
    title('3D Poisson equation');
end
end



function X=xx(A,pcase)
% Given A, a square SPD centrosymmetric matrix of size n, this algorithm produces
% double-cone matrices P, X and Y such that A=X*X' by using pchol (chol) factorization of
% diagonal blocks in similarity transformation of the matrix.
% if pcase=0 it uses chol
% if case=1 it uses pchol
% Author: Sarah Nataj, email:sarah.nataj@gmail.com
n=size(A,1);
% finding the similarity transformation
p = floor(n/2);
q = ceil(n/2);
pq = p+1:q;
A1 = A(1:p,1:p);
A2 = A(1:p,q+1:end);
A2 = A2(:,end:-1:1);
u = A(1:p,pq);
v = A(pq,1:p);
a = A(pq,pq);
B1 = [ A1+A2 , u*sqrt(2) ; v*sqrt(2), a ];
B2 = A1-A2;
% finding Cholesky factorization of the diagonal blocks
if pcase==0
    U1=chol(B1);
    U2=chol(B2);
elseif pcase==1
    U1=pchol(B1);
    U2=pchol(B2);
end
L1=U1';L2=U2';
% making X such that  A=XX'
d=q-p;
l=repmat(L1(end,1:end-d),d,1)/sqrt(2);
kap=repmat(L1(end,end),d,d);
L1=L1(1:end-d,1:end-d);
Lp=(L1+L2)/2;
Lm=(L1-L2)/2;
Zl=zeros(size(L1,2),d);
X=[ Lp, Zl, Lm(:,end:-1:1); l,kap, l(end:-1:1);  Lm(end:-1:1,:), Zl, Lp(end:-1:1,end:-1:1)];
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
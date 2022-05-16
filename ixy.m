function [X,Y]=ixy(A)
% Given A, a square centrosymmetric matrix, this algorithm produces two
% double-cone matrices X and Y such that A~X*Y by using ilu factorization of
% diagonal blocks in similarity transformation of the matrix.
% Author: Sarah Nataj, email:sarah.nataj@gmail.com
n=size(A,1);
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
% finding iLU factorization of the diagonal blocks
[L1,U1]=ilu(B1);
[L2,U2]=ilu(B2);
% making X and Y, such that A=~XY
d=q-p;
l=repmat(L1(end,1:end-d),d,1)/sqrt(2);
u=repmat(U1(1:end-d,end),1,d )/sqrt(2);
kap=repmat(L1(end,end),d,d);
rho=repmat(U1(end,end),d,d);
U1=U1(1:end-d,1:end-d);
L1=L1(1:end-d,1:end-d);
Up=(U1+U2)/2;
Um=(U1-U2)/2;
Lp=(L1+L2)/2;
Lm=(L1-L2)/2;
Zu=zeros(d,size(L1,2));
Zl=zeros(size(L1,2),d);
X=[ Lp, Zl, Lm(:,end:-1:1); l,kap, l(end:-1:1);  Lm(end:-1:1,:), Zl, Lp(end:-1:1,end:-1:1)];
Y=[ Up,  u, Um(:,end:-1:1); Zu,rho, Zu; Um(end:-1:1,:), u(end:-1:1), Up(end:-1:1,end:-1:1)];
end
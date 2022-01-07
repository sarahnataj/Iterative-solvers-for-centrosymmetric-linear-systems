function [P,X,Y]=xy(A,pcase)
% Given A, an  square centrosymmetric matrix of size n, this algorithm produces
% double-cone matrices P, X and Y such that A=P'*X*Y by using plu (lu) factorization of
% diagonal blocks in similarity transformation of the matrix.
% if pcase=0 it uses lu
% if pcase=1 it uses plu
% Author: Sarah Nataj, email:sarah.nataj@gmail.com
n=size(A,1);
%finding the similarity transformation
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
% finding LU factorization of the diagonal blocks
if pcase==0
    [L1,U1,P1]=lu(B1);
    [L2,U2,P2]=lu(B2);
elseif pcase==1
    [L1,U1,P1]=plu(B1);
    [L2,U2,P2]=plu(B2);
end
% making P, X and Y such that A=P'XY
d=q-p;
t=repmat(P1(end,1:end-d),d,1)/sqrt(2);
s=repmat(P1(1:end-d,end),1,d )/sqrt(2);
l=repmat(L1(end,1:end-d),d,1)/sqrt(2);
u=repmat(U1(1:end-d,end),1,d )/sqrt(2);
gam=repmat(P1(end,end),d,d);
kap=repmat(L1(end,end),d,d);
rho=repmat(U1(end,end),d,d);
U1=U1(1:end-d,1:end-d);
L1=L1(1:end-d,1:end-d);
P1=P1(1:end-d,1:end-d);
Up=(U1+U2)/2;
Um=(U1-U2)/2;
Lp=(L1+L2)/2;
Lm=(L1-L2)/2;
Pp=(P1+P2)/2;
Pm=(P1-P2)/2;
Zu=zeros(d,size(L1,2));
Zl=zeros(size(L1,2),d);
X=[ Lp, Zl, Lm(:,end:-1:1); l,kap, l(end:-1:1);  Lm(end:-1:1,:), Zl, Lp(end:-1:1,end:-1:1)];
Y=[ Up,  u, Um(:,end:-1:1); Zu,rho, Zu; Um(end:-1:1,:), u(end:-1:1), Up(end:-1:1,end:-1:1)];
P=[ Pp,  s, Pm(:,end:-1:1); t,gam, t(end:-1:1); Pm(end:-1:1,:), s(end:-1:1), Pp(end:-1:1,end:-1:1)];
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
function [X, Y, P]=ixytp(A, drop)
% Given A, a square centrosymmetric matrix, this algorithm produces
% double-cone matrices P, X and Y such that PA~X*Y by using ilutp factorization of
% diagonal blocks in similarity transformation of the matrix.
% drop is drop tolerance in ilutp.
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
% finding ilutp factorization of the diagonal blocks
[L1, U1, P1]=ilu(sparse(B1),struct('type','ilutp','droptol',drop));
[L2, U2, P2]=ilu(sparse(B2),struct('type','ilutp','droptol',drop));
% making P, X and Y, such that PA=~XY
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


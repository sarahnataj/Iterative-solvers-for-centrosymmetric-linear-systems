function X=ixxt(A,drop)
% Given A, a square SPD centrosymmetric matrix, this algorithm produces a
% double-cone matrix X such that A~XX' by using ichol factorization with drop tolerance of
% diagonal blocks in similarity transformation of the matrix.
% drop is drop tolerance in ichol.
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
% finding Cholesky factorization of the diagonal blocks
L1=ichol(B1,struct('shape','lower','type','ict','droptol',drop,'michol','on'));
L2=ichol(B2,struct('shape','lower','type','ict','droptol',drop,'michol','on'));
% making X such that A~XX'
d=q-p;
l=repmat(L1(end,1:end-d),d,1)/sqrt(2);
kap=repmat(L1(end,end),d,d);
L1=L1(1:end-d,1:end-d);
Lp=(L1+L2)/2;
Lm=(L1-L2)/2;
Zl=zeros(size(L1,2),d);
X=[ Lp, Zl, Lm(:,end:-1:1); l,kap, l(end:-1:1);  Lm(end:-1:1,:), Zl, Lp(end:-1:1,end:-1:1)];
end
function [error, A]=error_PoissonRobin2D(n)
% Spectral Legendre collocation method for
% -\laplacian u=f(x,y),
% with Robin boundary condition
% du\dn+a(x,y)u=0
% Input: n is the number of interior collocation nodes.
%     
% Output: error=\|u-u_ex\|_{\infty}, uex is the exact solution of the 2D 
%         Poisson equation at the collocation points, u is the
%         approximated solution.
%         A is the second order spectral differentiation matrix.
% Author: Sarah Nataj, email:sarah.nataj@gmail.com
n1=n+1; D=zeros(n1);
j=1:n-2; d=.5*sqrt(j.*(j+2)./(j+.5)./(j+1.5));
node=eig(diag(d,1)+diag(d,-1)); % Legendre Lobatto pts
node=[1; flipud(node); -1];
legn=legendre(n,node); legn=legn(1,:); % L_n(nodes)
for i=1:n1 % calculate deriv matrix
    legni=legn(i); nodei=node(i);
    for j=1:n1
        if i ~= j D(i,j)=legni/legn(j)/(nodei-node(j)); end; end;
end
D(1,1)=n*n1/4; D(n1,n1)=-n*n1/4;
w=2/n/(n+1)./legn.^2; W=spdiags(w',0,n1,n1); 
d=zeros(n1^2,1);
a=@(x,y) 2+x;
d(1:n1)=a(1,node);
d(end-n:end)=a(-1,node);
for j=1:n-1
    d(j*n1+1)=a(node(j+1),1); d(j*n1+n1)=a(node(j+1),n1); end
AA=D'*W*D;C=spdiags(d,0,n1^2,n1^2);
A=kron(W,AA)+kron(AA,W)+C;
[X,Y]=meshgrid(node);
uex=exp(Y).*((1-Y.^2).^2).*((1-X.^2).^2); uex=reshape(uex',n1^2,1);
f=exp(Y).*((4-12*X.^2).*(1-Y.^2).^2+(3+8*Y-10*Y.^2-8*Y.^3-Y.^4).*(1-X.^2).^2);
f=reshape(f',n1^2,1);
f=kron(W,W)*f;
u=A\f;
error=norm(u-uex,inf);
%figure(1), mesh(X,Y,reshape(u,n1,n1)); figure(2), mesh(X,Y,reshape(uex,n1,n1));

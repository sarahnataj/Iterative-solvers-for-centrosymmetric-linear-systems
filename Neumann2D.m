function [A, uex]=Neumann2D(n)
% Spectral Legendre collocation method for
% -\laplacian u+u=f(x,y),
% with Neumann condition
% du\dn=0,
% Input: n is the number of interior collocation nodes,
% Output: A is the 2D second spectral differentiation matrix
%        uex is the exact solution of the problem
% Author: Sarah Nataj, email:sarah.nataj@gmail.com
n1=n+1; D=zeros(n1);
j=1:n-2; d=.5*sqrt(j.*(j+2)./(j+.5)./(j+1.5));
node=eig(diag(d,1)+diag(d,-1)); % Legendre Lobatto pts
node=[1; flipud(node); -1];
legn=legendre(n,node); legn=legn(1,:); % L_n(nodes)
for i=1:n1 % calculate deriv matrix
    legni=legn(i); nodei=node(i);
    for j=1:n1
        if i ~= j
            D(i,j)=legni/legn(j)/(nodei-node(j));
        end
    end
end
D(1,1)=n*n1/4; D(n1,n1)=-n*n1/4;
w=2/n/(n+1)./legn.^2;
W=spdiags(w',0,n1,n1); 
B=D'*W*D;
A=kron(B,W)+kron(W,B)+kron(W,W);
[X,Y]=meshgrid(node);
uex=cos(pi*Y).*((1-X.^2).^2);uex=reshape(uex,(n+1)^2,1);
% f=cos(pi*Y).*(4*(1-3*X.^2)+(pi^2+1)*(X.^2-1).^2);
% f=reshape(f,(n+1)^2,1);
% f=kron(W,W)*f;
% u=A\f;
% error=norm(u-uex);
% spy(A)
end
function A=HelmholtzRobin2D(n,k)
% Spectral Legendre collocation method for 2D Helmholtz equation
% -\laplacian u+k^2u=f(x,y),
% with Robin boundary condition
% Input: n is the number of interior collocation nodes,
%        k is the wave number
% Output: A is the associated spectral differentiation matrix
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
[X,Y]=meshgrid(node,node);
C=afun(X,Y);
A=kron(B,W)+kron(W,B)-(k^2)*kron(W,W)+C;
%F=kron(B,W)+kron(W,B)+C;G=kron(W,W);
%s=sort(eig(F,G)); s(400)
%plot(s,'*r')
end


function C=afun(X,Y)
m=size(X,1);
C=2+X+Y.^2;
C(2:m-1,2:m-1)=zeros(m-2,m-2);
C=reshape(C,m^2,1);
C=spdiags(C,0,m^2,m^2);
end

function [A, C]=PoissonRobin2D(n)
% Spectral Legendre collocation method for
% -\Delta u=f(x,y),
% with Robin condition
% du\dn+a(x,y)u=0, 
% a(x,y)=2+x or a(x,y)=1+e+x+y^2 
% Input: n is the number of interior collocation nodes,
% Output: A is the associated spectral differentiation matrix
%         C is the preconditioner  for A, solving the same problem with 
%         a= average of a(x,y) on boundary
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
[X,Y]=meshgrid(node,node);
C=afun(X,Y);
B=D'*W*D;
A=kron(W,B)+kron(B,W)+C;
C=pfun(n1);
C=kron(W,B)+kron(B,W)+C;
end



function M=afun(X,Y)
m=size(X,1);
M=1+0.00001+X+(Y.^2);
M(2:m-1,2:m-1)=zeros(m-2,m-2);
M=reshape(M,m^2,1);
M=spdiags(M,0,m^2,m^2);
end
function M=pfun(m)
a=5/3+0.00001;
M=ones(m,m);
M(2:m-1,2:m-1)=zeros(m-2,m-2);
M=a*M;
M=reshape(M,m^2,1);
M=spdiags(M,0,m^2,m^2); 
end

% function M=afun(X,Y)
% m=size(X,1);
% M=2+X;
% M(2:m-1,2:m-1)=zeros(m-2,m-2);
% M=reshape(M,m^2,1);
% M=spdiags(M,0,m^2,m^2); 
% end
% function M=pfun(m)
% a=2;
% M=ones(m,m);
% M(2:m-1,2:m-1)=zeros(m-2,m-2);
% M=a*M;
% M=reshape(M,m^2,1);
% M=spdiags(M,0,m^2,m^2); 
% end


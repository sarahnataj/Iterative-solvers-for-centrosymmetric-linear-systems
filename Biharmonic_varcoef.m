function [A, C, M, uex]=Biharmonic_varcoef(n,k)
% Calculate spectral differentiation (pseudospectral) matrix associated with
% 2D variable-cofficient biharmonic equation given by
% \Delta(a \Delta u)=f u\in[-1,1], with homogeneous Dirichlet boundary conditions 
% a(x,y)=1+kx^2y^2 for k1>=a(x,y)>=k0>0, k>=0.
% Input: n is the number of interior nodes,
%        k is the constant in a(x,y)
% Output: A is the associated 2D spectral differentition matrix 
%        C is the sparse part of A.
%        M is the C plus dense part of A approximated by finite difference
%        at the colocation points
% Author: Sarah Nataj, email:sarah.nataj@gmail.com
method=2; % 1=Legendre  2=Cheb
[D, node]=PSDirv(n,method);
D2=D*D; 
B=(diag(1-node.^2)*D^2-4*diag(node)*D-2*eye(n+1))*diag([0; 1./(1-node(2:n).^2); 0]);
[X,Y]=meshgrid(node,node);
S=1+k*((X.^2).*(Y.^2)); S=diag(reshape(S,(n+1)^2,1));S=sparse(S);
A=kron(D2,eye(n+1))*S*kron(B,eye(n+1))...
    +kron(eye(n+1),D2)*S*kron(eye(n+1),B)...
    +kron(eye(n+1),D2)*S*kron(B,eye(n+1))...
    +kron(D2,eye(n+1))*S*kron(eye(n+1),B);
in=[];
    for jj=1:n-1
        in=[in, jj*(n+1)+2:jj*(n+1)+n];
    end
A=A(in,in);
% calculate C, the sparse part of A
C=kron(D2,speye(n+1))*S*kron(B,speye(n+1))...
    +kron(speye(n+1),D2)*S*kron(speye(n+1),B);
C=C(in,in);
% calculate dense part of A approximated by finite difference at the colocation points
N=fdpoisson(node);
S=S(in,in);
M=C+kron(speye(n-1),N)*S*kron(N,speye(n-1))...
    +kron(N,speye(n-1))*S*kron(speye(n-1),N);

% calcualte uex
[X,Y]=meshgrid(node(2:n));
uex=(sin(pi*X).^2).*(sin(pi*Y).^2);
uex=reshape(uex',(n-1)^2,1);

end

function M=fdpoisson(node)
n=size(node,1)-1;
for j=1:size(node)-1
    h(j)=node(j)-node(j+1);
end
B=zeros(n-1,n-1);
B(1,1)=2/(h(2)*h(1));
B(1,2)=-2/(h(2)*(h(2)+h(1)));
for j=2:n-2
    B(j,j)=2/(h(j)*h(j+1));
    B(j,j+1)=-2/(h(j+1)*(h(j)+h(j+1)));
    B(j,j-1)=-2/(h(j)*(h(j)+h(j+1)));
end
B(n-1,n-2)=-2/(h(n-1)*(h(n-1)+h(n)));
B(n-1,n-1)=2/(h(n)*h(n-1));
M=B;
%M=kron(speye(n-1),B)+kron(B,speye(n-1));
end
function [A,uex]=PseudoSpectral2D(n,order)
% Calculate 2D second/fourth order spectral differentiation (pseudospectral) matrices.
% Input: n is the number of interior collocation nodes,
%        order=2 for second order spectral differentiation matrix,
%        order=4 for fourth order spectral differentiation matrix.
% Output: A is the second/fourth spectral differentiation matrices,
%         uex is the exact solution of the 2D Poisson/biharmonic equations
%         at the collocation points
% Author: Sarah Nataj, email:sarah.nataj@gmail.com
method=2; % 1=Legendre  2=Cheb
[D, node]=PSDirv(n,method);
[X,Y]=meshgrid(node(2:n), node(2:n));
if order==2
    %SECOND ORDER
    DD=D*D; A=kron(speye(n-1),DD(2:n,2:n))+kron(DD(2:n,2:n),speye(n-1));
    %the exact solution and f for second order
    uex=sin(pi*X).*sin(2*pi*Y); uex=reshape(uex',(n-1)^2,1);
%     f=-(5*pi*pi)*sin(pi*X).*sin(2*pi*Y);f=reshape(f',(n-1)^2,1);
%     error=norm(uex-A\f);
elseif order==4
    %FOURTH ORDER
    B=(diag(1-node.^2)*D^2-8*diag(node)*D-12*eye(n+1))*D^2*diag([0; 1./(1-node(2:n).^2); 0]); %4th deriv matrix 
    B=B(2:n,2:n); 
    E=(diag(1-node.^2)*D^2-4*diag(node)*D-2*eye(n+1))*diag([0; 1./(1-node(2:n).^2); 0]);%2nd deriv matrix
    E=E(2:n,2:n); 
    A=kron(B,speye(n-1))+kron(speye(n-1),B)+2*kron(E,E);%kron(E,I)*kron(I*E); or instead of E, use D2
    uex=(1+cos(pi*X)).*(1-2*Y.^2+Y.^4); uex=reshape(uex',(n-1)^2,1);   
%     f=pi^4*cos(pi*X).*(1-2*Y.^2+Y.^4)+24*(1+cos(pi*X))-8*pi^2*cos(pi*X).*(3*Y.^2-1); f=reshape(f',(n-1)^2,1);
%     error=norm(uex-A\f);
end
end
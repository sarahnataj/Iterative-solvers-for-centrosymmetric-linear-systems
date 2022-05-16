function [error, A]=error_singular_perturbation(n,e)
% spectral differentiation (pseudospectral) matrices
% singular perturbation problem
% ep*\Delta u+u_x=0.
% Input: n is the number of interior collocation nodes,
%        ep: epsilon
% Output: error=\|u-u_ex\|_{\infty}, uex is the exact solution of the
%         equation at the collocation points, u is the approximated solution.
%         A is the associated spectral differentiation matrix.
%         
% Author: Sarah Nataj, email:sarah.nataj@gmail.com
e2=e^2;
method=2; % 1=Legendre  2=Cheb
[D, node]=PSDirv(n,method);
[X,Y]=meshgrid(node(2:n), node(2:n));
%SECOND ORDER
DD=D*D;
A=-e2*kron(speye(n-1),DD(2:n,2:n))...
    -e2*kron(DD(2:n,2:n),speye(n-1))+kron(speye(n-1),D(2:n,2:n));
%the exact solution and f for second order
uex=(X+1).*(1-exp((X-1)/e)).*cos(pi*Y/2);
uex=reshape(uex',(n-1)^2,1);
f=e^2*((2*exp((X - 1)/e).*(cos((pi*Y)/2))/e)+...
    (exp((X - 1)/e).*cos((pi*Y)/2).*((X + 1))/e^2) -...
    (pi^2*cos((pi*Y)/2).*((exp((X - 1)/e) - 1)).*((X + 1))/4)) - ...
    cos((pi*Y)/2).*((exp((X - 1)/e) - 1)) - ...
    (exp((X - 1)/e).*(cos((pi*Y)/2)).*((X + 1))/e);
f=reshape(f',(n-1)^2,1);
u=A\f;
error=norm(uex-u,inf);
%figure (1); mesh(X,Y,reshape(uex,n-1,n-1))
%figure (2); mesh(X,Y,reshape(u,n-1,n-1))
%end
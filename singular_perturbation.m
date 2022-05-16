function [A,uex]=singular_perturbation(n,e)
% spectral differentiation (pseudospectral) matrix for 
% singular pertubation problem homogeneous Dirichlet boundary conditions
% ep*\Delta u+u_x=0, 
% Input: n is the number of interior collocation nodes,
%        ep: epsilon
% Output: A is the associated  spectral differentiation matrix,
%         uex is the exact solution of the equation at the collocation
%         points.
% Author: Sarah Nataj, email:sarah.nataj@gmail.com
e2=e^2;
method=2; % 1=Legendre  2=Cheb
[D, node]=PSDirv(n,method);
[X,Y]=meshgrid(node(2:n), node(2:n));
DD=D*D;
A=-e2*kron(speye(n-1),DD(2:n,2:n))...
    -e2*kron(DD(2:n,2:n),speye(n-1))+kron(speye(n-1),D(2:n,2:n));
%the exact solution and f for second order
uex=(X+1).*(1-exp((X-1)/e)).*cos(pi*Y/2);
uex=reshape(uex',(n-1)^2,1);

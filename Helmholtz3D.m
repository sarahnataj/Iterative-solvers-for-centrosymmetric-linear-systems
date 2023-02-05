function [A, f]=Helmholtz3D(n,a)
% Calculate 3D second order spectral differentiation (pseudospectral)
% matrix for 3D Helmholtz equation, -(\Delta +a)u=f with with homogeneous Dirichlet boundary conditions
% Input: n is the number of interior collocation nodes,
%        a is the square of the wave number.
% Output: A is the associated spectral differentiation matrix for Helmholtz
%        equation.
%        uex is the exact solution of 3D Helmholtz equation at
%        the collocation points
% Author: Sarah Nataj, email:sarah.nataj@gmail.com
method=2; % 1=Legendre  2=Cheb
[D, node]=PSDirv(n,method);
[X,Y,Z]=meshgrid(node(2:n), node(2:n),node(2:n));
DD=D*D;
A=-kron(kron(speye(n-1),DD(2:n,2:n)),speye(n-1))- ...
    kron(kron(DD(2:n,2:n),speye(n-1)),speye(n-1))-...
    kron(kron(speye(n-1),speye(n-1)),DD(2:n,2:n))-...
    (a*speye((n-1)^3));
uex=sin(pi*X).*sin(pi*Y).*sin(pi*Z); uex=reshape(uex,(n-1)^3,1);
f=(3*pi^2-a)*uex;
end
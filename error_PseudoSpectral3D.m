function [error, A]=error_PseudoSpectral3D(n)
% Calculate second order spectral differentiation (pseudospectral) 
% matrix for 3D Poisson equation with homogeneous Dirichlet boundary conditions.
% Input: n is the number of interior collocation nodes.
%     
% Output: error=\|u-u_ex\|_{\infty}, uex is the exact solution of the 3D 
%         Poisson equation at the collocation points, u is the
%         approximated solution.
%         A is the second order spectral differentiation matrix.
% Author: Sarah Nataj, email:sarah.nataj@gmail.com
method=2; % 1=Legendre  2=Cheb
%n=21;
[D, node]=PSDirv(n,method);
DD=D*D; A=kron(kron(speye(n-1),DD(2:n,2:n)),speye(n-1))+ ...
    kron(kron(DD(2:n,2:n),speye(n-1)),speye(n-1))+...
    kron(kron(speye(n-1),speye(n-1)),DD(2:n,2:n));
[X,Y,Z]=meshgrid(node(2:n), node(2:n),node(2:n));
uex=sin(pi*X).*sin(pi*Y).*sin(pi*Z); uex=reshape(uex,(n-1)^3,1);
f=-3*pi^2*uex;
error=norm(uex-A\f,inf);
end

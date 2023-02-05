function [error, A]=error_Poisson_varcoef(n,dim,k, acase)
% Calculate 2D spectral differentiation (pseudospectral) matrix associated with
% variable-cofficient Poisson operator,
% -div(a grad u)=f u\in[-1,1], u=0 in boundary, a(x,y)=1+kx^2y^2 for
% k1>=a(x,y)>=k0>0, k>=0
% Input: n is the number of interior collocation nodes
%        dim is the dimension,
%        dim=1 for 1D problem -(a u_x)_x=f, a=1+kx^2,  u(-1)=0=u(1)
%        dim=2 for 2D problem -div(a grad u)=f 
%        k is the constant in a(x,y)=1+kx^2y^2 
%        acase is to choose a(x,y)
% Output: error=\|u-u_ex\|_{\infty}, uex is the exact solution of the 1D/2D 
%         diffusion equations at the collocation points, u is the
%         approximated solution.
%         A is the associated spectral differentiation matrix.
% Author: Sarah Nataj, email:sarah.nataj@gmail.com
method=2; % 1=Legendre  2=Cheb
[D, node]=PSDirv(n,method);
if dim==1 
    X=node(1:end);
    S=1+k*(X.^2);S=diag(S);
    D2=D*S*D;
    A=-D2(2:n,2:n);
    X=node(2:n);
    uex=sin(pi*X);
    f=pi^2*(1+k*(X.^2)).*(sin(pi*X))-2*pi*k*X.*cos(pi*X);
    u=A\f;
    norm(uex-u,inf)
elseif dim==2 
    [X,Y]=meshgrid(node,node);
    if acase==1
    S=1+k*((X.^2).*(Y.^2)); 
    else
       S=1+k*((X.^2).*(Y.^4));
    end
    S=diag(reshape(S,(n+1)^2,1));S=sparse(S);
    A=-kron(D,speye(n+1))*S*kron(D,speye(n+1))-kron(speye(n+1),D)*S*kron(speye(n+1),D);
    in=[];
    for jj=1:n-1
        in=[in, jj*(n+1)+2:jj*(n+1)+n];
    end
    A=A(in,in);
    [X,Y]=meshgrid(node(2:n));
    uex=sin(pi*X).*sin(pi*Y);uex=reshape(uex',(n-1)^2,1);
    if acase==1
    f=2*(1+k*((X.^2).*(Y.^2))).*(pi*pi).*sin(pi*X).*sin(pi*Y)-...
            2*k*pi*(X.*Y).*(Y.*cos(pi*X).*sin(pi*Y)+X.*sin(pi*X).*cos(pi*Y));
    else
f=-2*pi*k*X.*(Y.^4).*cos(pi*X).*sin(pi*Y)...
    -4*pi*k*(X.^2).*(Y.^3).*sin(pi*X).*cos(pi*Y)...
+2*(pi^2)*(1+k*((X.^2).*(Y.^4))).*sin(pi*X).*sin(pi*Y);
    end
    f=reshape(f,(n-1)^2,1);
    u=A\f;
    error=norm(uex-u,inf);
end

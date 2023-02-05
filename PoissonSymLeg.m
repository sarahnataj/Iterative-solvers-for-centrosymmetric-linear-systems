function [A,uex]=PoissonSymLeg(n)
% Calculate 2D second order symmetric positive definite (SPD) spectral
% differentiation (pseudospectral) matrix using Legendre Gauss Lobatto points.
%
% Input: n is the number of interior collocation nodes
% Output: A is the associated SPD spectral differentiation matrix
%        uex is the exact solution of the Poisson equation at the
%        collocation points.
% Author: Sarah Nataj, email:sarah.nataj@gmail.com
n1=n+1; D=zeros(n1);
j=1:n-2; d=.5*sqrt(j.*(j+2)./(j+.5)./(j+1.5));
node=eig(diag(d,1)+diag(d,-1));
node=[1; flipud(node); -1];
legn=legendre(n,node); legn=legn(1,:);
for i=1:n1
    legni=legn(i); nodei=node(i);
    for j=1:n1
        if i ~= j
            D(i,j)=legni/legn(j)/(nodei-node(j));
        end
    end
end
D(1,1)=n*n1/4; D(n1,n1)=-n*n1/4;
[X,Y]=meshgrid(node(2:n));
w=2/n/(n+1)./legn.^2;
C=diag(w)*D; C=C'*diag(1./w)*C; C=C(2:n,2:n);
ww=sqrt(w(2:n))'; M=diag(1./ww)*C*diag(1./ww);
A=kron(M,speye(n-1))+kron(speye(n-1),M);
uex=sin(pi*X).*sin(2*pi*Y); uex=reshape(uex',(n-1)^2,1);
%  f=(5*pi*pi)*sin(pi*X).*sin(2*pi*Y);f=reshape(f',(n-1)^2,1);
%  f=kron(diag(ww),diag(ww))*f;
%  vh=A\f;
%  uh=kron(diag(1./ww),diag(1./ww))*vh;
%  error=norm(uex-uh);
end


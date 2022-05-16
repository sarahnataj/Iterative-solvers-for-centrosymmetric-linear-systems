function [Nnop,Nilu,Nixy,Nilutp,Nixytp]=gmres_poisson(n,ecase)
% GMRES preconditioned by incomplete double-cone factorization
% for 2D/3D Poisson equation and
%     2D Poisson equation with variable coefficients (diffusion equation).
%     2D singular perturbation problem
% Input: n is the number of the collocation nodes
%       toltp = ilutp drop tolerance
%       tol = GMRES tolerance
% Output: The number of iterations for GMRES and density of proposed preconditioners.
% Author: Sarah Nataj, email:sarah.nataj@gmail.com
tol = 1e-6;
toltp=1e-3;
if ecase==1
A=PseudoSpectral2D(n,2);
elseif ecase==2
A=PseudoSpectral3D(n);
elseif ecase==3
A=Poisson_varcoef(n,2,10);
elseif ecase==4
[A,uex]=singular_perturbation(n,0.01);
J=flipud(speye(size(A,1)));
A=0.5*(A+J*A*J);
%A=0.5*(A+A');
end
x0 = zeros(size(A,1),1);
b=randn(size(A,1),1);
maxit = size(A,1);
[x0,fl0,rr0,it0,rv0] = gmres(A,b,[],tol,maxit,[],[],x0);
%size(full(A),1)
Nnop=it0(2);
%
[L1, U1]=ilu(A);
[x1,fl1,rr1,it1,rv1] = gmres(A,b,[],tol,maxit,L1,U1);
Nilu=it1(2);
%
[X1, Y1]=ixy(A);
[x2,fl2,rr2,it2,rv2]= gmres(A,b,[],tol,maxit,X1,Y1);
Nixy=it2(2);
%
[L2, U2,P]=ilu(A,struct('type','ilutp','droptol',toltp));
[x3,fl3,rr3,it3,rv3]= gmres(P*A,P*b,[],tol,maxit,L2,U2);
Nilutp=it3(2);
%
[X2, Y2, P1]=ixytp(A,toltp);
[x4,fl4,rr4,it4,rv4]= gmres(P1*A,P1*b,[],tol,maxit,X2,Y2);
Nixytp=it4(2);
end




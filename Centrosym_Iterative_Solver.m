% This program applys GMRES/PCG preconditioned by incomplete double-cone 
% factorization for centrosymetric matrices arising from spectral 
% discretization of 
%     2D/3D Poisson equations with homogeneous Dirichlet and Robin boundary conditions.
%     2D diffusion equation with homogeneous Dirichlet boundary conditions.
%     2D Nuamann problem -\laplacian u+u=f(x,y) with Neumann boundary conditions.
%     2D biharmonic equation with homogeneous Dirichlet boundary conditions.
%     2D variable-coefficients biharmonic equation with homogeneous Dirichlet boundary conditions.
%     2D/3D Helmholtz equations with homogeneous Dirichlet and Robin boundary conditions.
%
% Input: n is the number of the collocation nodes
%       toltp = ilutp drop tolerance
%       tol = GMRES tolerance
%       ecase= the equation
%       precase= the preconditioner 
%
% The preconditioners are discussed in the manuscript: Structure-preserving 
% solvers for centrosymmetric linear systems with applications to spectral 
% methods, by Chen Greif, Sarah Nataj and Manfred Trummer.
%
% Author: Sarah Nataj, email:sarah.nataj@gmail.com
clc;
clear all;
close all;
tol = 1e-6;
toltp=1e-3;
n=19;
ecase=11;
if ecase==1
    disp('2D Poisson equation with homogeneous Dirichlet boundary conditions');
    Precond_Poisson(n,1,tol,toltp)
elseif ecase==2
    disp('3D Poisson equation with homogeneous Dirichlet boundary conditions');
    Precond_Poisson(n,2,tol,toltp)
elseif ecase==3
    disp('2D diffusion equation with homogeneous Dirichlet boundary conditions')
    Precond_Poisson(n,3,tol,toltp)
elseif ecase==4
    disp('Symmetric Legendre spectral matrix for 2D Poisson equation with homogeneous Dirichlet boundary conditions')
    Precond_Possion_SPD(n,0,tol,toltp)
elseif ecase==5
    disp('Symmetric Legendre spectral matrix for 2D Neumann problem')
    Precond_Possion_SPD(n,1,tol,toltp)
elseif  ecase==6
    disp('Finite difference for 2D Poisson equation with homogeneous Dirichlet boundary conditions')
    Precond_Possion_SPD(n,2,tol,toltp)
elseif  ecase==7
    disp('Finite difference for 3D Poisson equation with homogeneous Dirichlet boundary conditions')
    Precond_Possion_SPD(n,3,tol,toltp)
elseif  ecase==8
    disp('2D Poisson equation with Robin bundary conditions');
    precase=1;% 1=centrosymmetric part, 2=approximated matrix
    Precond_Possion_Robin(n,precase,tol,toltp)
elseif ecase==9
    disp('3D Helmholtz equation with homogeneous Dirichlet boundary conditions')
    a=10;
    Precond_Helmholtz3D(n,a,tol,toltp)
elseif ecase==10
    disp('2D Helmholtz equation with Robin boundary conditions')
    a=sqrt(8.4255386e+03);
    Precond_HelmholtzRobin2D(n,a, tol, toltp)
elseif ecase==11
    disp('2D biharmonic equation with homogeneous Dirichlet boundary conditions');
    precase=3;% 1=G, 2=C, 3=M
    Precond_Biharmonic(n,0, precase, tol,toltp)
    elseif ecase==12
    disp('2D variable coeficients biharmonic equation with homogeneous Dirichlet boundary conditions');
    k=10;
    precase=3;%  2=C, 3=M
    Precond_Biharmonic(n,k, precase, tol,toltp)
end
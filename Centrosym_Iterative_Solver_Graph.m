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
% Input:toltp = ilutp drop tolerance
%       tol = GMRES tolerance
%       ecase= the equation
%
% We used GMRES with restart for 3D Helmholtz equation, for other PDEs, GMRES
% without restart have been used.
%
% The results have been proposed in the manuscript:
% Incomplete double-cone factorizations of centrosymmetric matrices arising
% in spectral methods by Chen Greif, Sarah Nataj and Manfred Trummer.
%
% Author: Sarah Nataj, email:sarah.nataj@gmail.com Feb 2023
clc;
clear all;
close all;
tol = 1e-6;
toltp=1e-4;
ecase=4;
if ecase==1
    disp('2D Poisson equation with homogeneous Dirichlet boundary conditions');
    Precond_Poisson(20,1,tol,toltp)
elseif ecase==2
    disp('3D Poisson equation with homogeneous Dirichlet boundary conditions');
    Precond_Poisson(18,2,tol,toltp)
elseif ecase==3
    disp('2D diffusion equation with homogeneous Dirichlet boundary conditions')
    Precond_Poisson(21,3,tol,toltp)
elseif ecase==4
    disp('Symmetric Legendre spectral matrix for 2D Poisson equation with homogeneous Dirichlet boundary conditions')
    Precond_Possion_SPD(21,0,tol,toltp)
elseif ecase==5
    disp('Symmetric Legendre spectral matrix for 2D Neumann problem')
    Precond_Possion_SPD(19,1,tol,toltp)
elseif  ecase==6
    disp('Finite difference for 2D Poisson equation with homogeneous Dirichlet boundary conditions')
    Precond_Possion_SPD(51,2,tol,toltp)
elseif  ecase==7
    disp('Finite difference for 3D Poisson equation with homogeneous Dirichlet boundary conditions')
    Precond_Possion_SPD(31,3,tol,toltp)
elseif  ecase==8
    disp('2D Poisson equation with Robin bundary conditions');
    precase=2;% 1=centrosymmetric part, 2=approximated matrix
    Precond_Possion_Robin(21,precase,tol,toltp)
elseif ecase==9
    disp('3D Helmholtz equation with homogeneous Dirichlet boundary conditions')
    %a=[219.5988, 605.2396,1052.7138, 2195.9064, 2434.9594,3362.1876,4096.6722, 9794.2739, 18759.0058];
    a=9794.27;
    %a=1000;
    Precond_Helmholtz3D(15,a,tol)
elseif ecase==10
    disp('2D Helmholtz equation with Robin boundary conditions')
    a=sqrt(8.4255386e+03);
    Precond_HelmholtzRobin2D(20,a, tol, toltp)
elseif ecase==11
    disp('2D biharmonic equation with homogeneous Dirichlet boundary conditions');
    precase=3;% 1=G, 2=C, 3=M
    Precond_Biharmonic(21,0, precase, tol,toltp)
elseif ecase==12
    disp('2D variable coeficients biharmonic equation with homogeneous Dirichlet boundary conditions');
    k=10;
    precase=2;%  2=C, 3=M
    Precond_Biharmonic(24,k, precase, tol,toltp)
end
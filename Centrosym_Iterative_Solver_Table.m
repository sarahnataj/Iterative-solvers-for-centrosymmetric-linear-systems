% This program results number of iterations versus n for GMRES preconditioned by incomplete double-cone
% factorization for centrosymetric matrices arising from spectral discretization of
%     2D/3D Poisson equations with homogeneous Dirichlet and Robin boundary conditions.
%     2D diffusion equation with homogeneous Dirichlet boundary conditions.    
%     2D biharmonic equation with homogeneous Dirichlet boundary conditions.
%     2D variable-coefficients biharmonic equation with homogeneous Dirichlet Dirichlet boundary conditions.
%
% Input: ecase= the equation
%
% The results have been proposed in the manuscript:
% Incomplete double-cone factorizations of centrosymmetric matrices arising
% in spectral methods by Chen Greif, Sarah Nataj and Manfred Trummer.
%
% Author: Sarah Nataj, email:sarah.nataj@gmail.com Feb 2023
clc
clear all;
close all;
ecase=7;%choose from 1 to 7

if ecase==1
    N=10:2:20;
    disp('2D Poisson equation with homogeneous Dirichlet boundary conditions');
elseif ecase==2
    N=8:2:16;
    disp('3D Poisson equation with homogeneous Dirichlet boundary conditions');
elseif ecase==3
    N=10:2:20;
    disp('2D diffusion equation with homogeneous Dirichlet boundary conditions')
elseif ecase==4
    N=8:2:18;
    disp('2D biharmonic equation with homogeneous Dirichlet boundary conditions');
    disp('Preconditioned by the sparse part of the biharmonic operator');
elseif ecase==5
    N=10:1:18;
    disp('2D biharmonic equation with homogeneous Dirichlet boundary conditions');
    disp(['Preconditioned by the sparse part + dense part of' ...
        ' A approximated by finite difference at the collocation points.']);
elseif ecase==6
    N=8:2:18;
    disp('2D variable-coefficient biharmonic equation with homogeneous Dirichlet boundary conditions');
    disp('Preconditioned by the sparse part of the biharmonic operator');
elseif ecase==7
    N=8:2:18;
    disp('2D variable-coefficient biharmonic equation with homogeneous Dirichlet boundary conditions');
    disp(['Preconditioned by the sparse part + dense part of' ...
        ' A approximated by finite difference at the collocation points.'])
end

for i=1:size(N,2)
    if ecase==1
        [N_nop(i),N_ilu(i),N_ixy(i),N_ilutp(i),N_ixytp(i)]=gmres_poisson(N(i),1);
    elseif ecase==2
        [N_nop(i),N_ilu(i),N_ixy(i),N_ilutp(i),N_ixytp(i)]=gmres_poisson(N(i),2);
    elseif ecase==3
        [N_nop(i),N_ilu(i),N_ixy(i),N_ilutp(i),N_ixytp(i)]=gmres_poisson(N(i),3);
    elseif ecase==4
        [N_nop(i),N_ilu(i),N_ixy(i),N_ilutp(i),N_ixytp(i)]=gmres_biharmonic(N(i),0,2);
    elseif ecase==5
        [N_nop(i),N_ilu(i),N_ixy(i),N_ilutp(i),N_ixytp(i)]=gmres_biharmonic(N(i),0,3);
    elseif ecase==6
        [N_nop(i),N_ilu(i),N_ixy(i),N_ilutp(i),N_ixytp(i)]=gmres_biharmonic(N(i),10,2);
    elseif ecase==7
        [N_nop(i),N_ilu(i),N_ixy(i),N_ilutp(i),N_ixytp(i)]=gmres_biharmonic(N(i),10,3);
    end
end
Results=table(N',N_nop',N_ilu',N_ixy',N_ilutp',N_ixytp',...
        'VariableNames',{'nodes', 'no precond','ILU(0)','IXY(0)', 'ILUTP','IXYTP'})

semilogy(N,N_nop,'-o',N,N_ilu,'-o',N,N_ixy,'-o',N,N_ilutp,'r-o',N,N_ixytp','b-o','LineWidth',2)
legend('N_nop','N_ilu','N_ixy','N_ilutp','N_ixytp','Location','best')
title('GMRES convergence without and with proposed preconditioners');
xlabel('N')
ylabel('Iteration number')
plotformat(1.5,6)
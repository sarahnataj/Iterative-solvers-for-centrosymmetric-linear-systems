% This program examines the accuracy of the spectral discretization methods
% which have been proposed in the manuscript:
% Structure-preserving solvers for centrosymmetric linear systems with
% applications to spectral methods, by Chen Greif, Sarah Nataj and Manfred
% Trummer.
%
% Input: ecase is the equation. ecase=1,2,.., 12
% Output: The figure for error \|u-u_ex\|_{\infty} versus N, the number of
%             interior collocation points. We use Matlab backslash operator 
%             A\b to calculate the approximate solution u. 
%         The number of interior collocation nodes, N, and size of the 
%             associated matrix where the error is the smallest.                
%         The sparsity of A: 1= it is sparse, 0= it is dense.
%         The centrosymmetry of A: 0= it is not centrosymmetric
%                                  1= it is centrosymmetric with error<10^-8
%                                  2= it is centrosymmetric with error<10^-4
%         condition number of A, where is the error is the smallest
%
% Author: Sarah Nataj, email:sarah.nataj@gmail.com
clear all
close all
clc
ecase=3;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ecase==1
    N=10:2:22;
    disp(['1D Poisson equation with homogeneous Dirichlet boundary conditions,' ...
        ' accuracy of Chebyshev collocation method']);
    for i=1:size(N,2)
        [e(i), A]=error_PseudoSpectral1D(N(i),2);
    end
    interior_nodes=N(i-1)
    error=e(i-1)
    sizeA=size(A)
    sparsity=issparse(A)
    centrosymmetry=iscentrosym(A)
    condition=cond(A)
elseif ecase==2
    N=10:2:28;
    disp(['2D Poisson equation with homogeneous Dirichlet boundary conditions' ...
        ', accuracy of Chebyshev collocation method']);
    for i=1:size(N,2)
        [e(i), A]=error_PseudoSpectral2D(N(i),2);
    end
    interior_nodes=N(i-1)
    error=e(i-1)
    sizeA=size(A)
    sparsity=issparse(A)
    centrosymmetry=iscentrosym(A)
    condition=cond(full(A))
elseif ecase==3
    N=10:2:22;
    disp(['3D Poisson equation with homogeneous Dirichlet boundary conditions,' ...
        ' accuracy of Chebyshev  collocation method']);
    for i=1:size(N,2)
        [e(i), A]=error_PseudoSpectral3D(N(i));
    end
   interior_nodes=N(i-1)
    error=e(i-1)
    sizeA=size(A)
    sparsity=issparse(A)
    centrosymmetry=iscentrosym(A)
    condition=cond(full(A))
elseif ecase==4
    N=10:2:22;
    disp(['2D diffusion equation with homogeneous Dirichlet boundary conditions,' ...
        ' accuracy of Chebyshev collocation method']);
    for i=1:size(N,2)
        [e(i), A]=error_Poisson_varcoef(N(i),2,10,1);
    end
    interior_nodes=N(i-1)
    error=e(i-1)
    sizeA=size(A)
    sparsity=issparse(A)
    centrosymmetry=iscentrosym(A)
    condition=cond(full(A))
elseif ecase==5
    N=8:2:18;
    disp(['1D biharmonic equation with homogeneous Dirichlet boundary conditions,' ...
        ' accuracy of Chebyshev collocation method']);
    for i=1:size(N,2)
        [e(i), A]=error_PseudoSpectral1D(N(i),4);
    end
    interior_nodes=N(i-1)
    error=e(i-1)
    sizeA=size(A)
    sparsity=issparse(A)
    centrosymmetry=iscentrosym(A)
    condition=cond(A)
elseif ecase==6
    N=10:2:24;
    disp(['2D biharmonic equation with homogeneous Dirichlet boundary conditions,' ...
        ' accuracy of Chebyshev collocation method']);
    for i=1:size(N,2)
        [e(i), A]=error_Biharmonic_varcoef(N(i),0,1);
    end
    interior_nodes=N(i-1)
    error=e(i-1)
    sizeA=size(A)
    sparsity=issparse(A)
    centrosymmetry=iscentrosym(A)
    condition=cond(A)
elseif ecase==7
    N=10:2:26;
    disp(['2D variable-coefficient biharmonic equation with homogeneous Dirichlet boundary conditions,' ...
        ' accuracy of Chebyshev collocation method']);
    for i=1:size(N,2)
        [e(i), A]=error_Biharmonic_varcoef(N(i),10,1);
    end
    interior_nodes=N(i-1)
    error=e(i-1)
    sizeA=size(A)
    sparsity=issparse(A)
    centrosymmetry=iscentrosym(A)
    condition=cond(A)
elseif ecase==8
    N=10:2:22;
    disp(['3D Helmholtz equation with homogeneous Dirichlet boundary conditions,' ...
        ' accuracy of Chebyshev collocation method']);
    for i=1:size(N,2)
        [e(i), A]=error_Helmholtz3D(N(i),100);
    end
    interior_nodes=N(i-1)
    error=e(i-1)
    sizeA=size(A)
    sparsity=issparse(A)
    centrosymmetry=iscentrosym(A)
    condition=cond(full(A))
elseif ecase==9
    N=10:2:26;
    disp(['2D Poisson equation with homogeneous Dirichlet boundary conditions' ...
        ', accuracy of symmetric Legendre collocation method']);
    for i=1:size(N,2)
        [e(i), A]=error_PoissonSymLeg(N(i));
    end
    interior_nodes=N(i-1)
    error=e(i-1)
    sizeA=size(A)
    sparsity=issparse(A)
    centrosymmetry=iscentrosym(A)
    condition=cond(full(A))
elseif ecase==10
    N=10:2:20;
    disp(['\Delta u+u=f(x,y) with Neumann boundary conditions,' ...
        ' accuracy of symmetric Legendre collocation method']);
    for i=1:size(N,2)
        [e(i), A]=error_Neumann2D(N(i));
    end
    interior_nodes=N(i-1)
    error=e(i-1)
    sizeA=size(A)
    sparsity=issparse(A)
    centrosymmetry=iscentrosym(A)
    condition=cond(full(A))
elseif ecase==11
    N=10:2:18;
    disp(['2D Poisson equation with Robin boundary conditions,' ...
        ' accuracy of symmetric Legendre collocation method']);
    for i=1:size(N,2)
        [e(i), A]=error_PoissonRobin2D(N(i));
    end
    interior_nodes=N(i-1)
    error=e(i-1)
    sizeA=size(A)
    sparsity=issparse(A)
    centrosymmetry=iscentrosym(A)
    condition=cond(full(A))
elseif ecase==12
    N=80:2:90;
    disp(['2D singular perturbation problem  with homogeneous Dirichlet boundary conditions,' ...
        ' accuracy of Chebyshev collocation method']);
    for i=1:size(N,2)
        [e(i), A]=error_singular_perturbation(N(i),0.01);
    end
    interior_nodes=N(i-1)
    error=e(i-1)
    sizeA=size(A)
    sparsity=issparse(A)
    centrosymmetry=iscentrosym(A)
    condition=cond(full(A))

end

semilogy(N,e,'-o','LineWidth',2)
if  ecase==1
    title('1D Poisson equation');
elseif ecase==2
    title('2D Poisson equation');
elseif ecase==3
    title('3D Poisson equation');
elseif ecase==4
    title('2D diffusion equation');
elseif ecase==5
    title('1D biharmonic equation');
elseif ecase==6
    title('2D biharmonic equation');
elseif ecase==7
    title('2D variable-coefficient biharmonic equation');
elseif ecase==8
    title('3D Helmholtz equation with Dirichlet boundary conditions');
elseif ecase==9
    title('2D Poisson equation with Dirichlet boundary condition');
elseif ecase==10
    title('\Delta u+u=f(x,y) with Neumann boundary conditions');
elseif ecase==11
    title('2D Poisson equation with Robin boundary condition');
elseif ecase==12
    title('2D singular perturbation problem');

end
xlabel('N')
ylabel('error')
plotformat(1.5,6)

function c=iscentrosym(A)
J=flipud(eye(size(A,1)));B=0.5*(A+J*A*J);
diffnorm=norm(A-B);
if diffnorm<=1e-8
    c=1;
elseif diffnorm<=1e-4
    c=2;
else
    c=0;
end
end
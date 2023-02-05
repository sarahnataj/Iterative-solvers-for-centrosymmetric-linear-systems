% This program examines the accuracy of the spectral discretization
% methods for common linear PDEs: 
%     2D/3D Poisson equations with homogeneous Dirichlet and Robin boundary conditions.
%     2D diffusion equation with homogeneous Dirichlet boundary conditions.
%     2D Nuamann problem -\laplacian u+u=f(x,y) with Neumann boundary conditions.
%     2D biharmonic equation with homogeneous Dirichlet boundary conditions.
%     2D variable-coefficients biharmonic equation with homogeneous Dirichlet boundary conditions.
%     3D Helmholtz equations with homogeneous Dirichlet and Robin boundary conditions.
%
% Input:  ecase= the equation. ecase=1,2,..
% Output: The figure for error \|u-u_ex\|_{\infty} versus N, the number of
%             collocation points. We use Matlab backslash operator 
%             A\b to calculate the approximate solution u. 
%         The number of collocation nodes, N
%         Error \|u-u_ex\|_{\infty}
%         Size of the associated square matrix                 
%         The sparsity of A: true = it is sparse, false= it is dense
%         The centrosymmetry of A: true= it is centrosymmetric
%                                  false= it is not centrosymmetric with error>10^-8
%         Condition number of A
%
% The results have been proposed in the manuscript:
% Incomplete double-cone factorizations of centrosymmetric matrices arising
% in spectral methods by Chen Greif, Sarah Nataj and Manfred Trummer.
%
% Author: Sarah Nataj, email:sarah.nataj@gmail.com, Feb 2023
clear all
close all
clc
ecase=9;%choose from 1 to 11
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ecase==1
    N=16:1:24;
    disp('Accuracy of Chebyshev collocation method for 1D Poisson equation with homogeneous Dirichlet boundary conditions');
    for i=1:size(N,2)
        [e(i), A]=error_PseudoSpectral1D(N(i),2);
        sizeA(i)=size(A,1);
        sparsity(i)=issparse(A);
        centrosymmetry(i)=iscentrosym(A);
        condition(i)=cond(full(A));
    end
    Results=table(N',e',sizeA',sparsity',centrosymmetry',condition',...
        'VariableNames',{'number of nodes', 'error',....
        'size of the matrix','sparsity', 'centrosymmetry','condition number'})
elseif ecase==2
    N=18:2:28;
    disp('Accuracy of Chebyshev collocation method for 2D Poisson equation with homogeneous Dirichlet boundary conditions');
    for i=1:size(N,2)
        [e(i), A]=error_PseudoSpectral2D(N(i)-1,2);
        sizeA(i)=size(A,1);
        sparsity(i)=issparse(A);
        centrosymmetry(i)=iscentrosym(A);
        condition(i)=cond(full(A));
    end
    Results=table(N',e',sizeA',sparsity',centrosymmetry',condition',...
        'VariableNames',{'number of nodes', 'error',....
        'size of the matrix','sparsity', 'centrosymmetry','condition number'})
elseif ecase==3
    N=21:1:22;
    disp('Accuracy of Chebyshev collocation method for 3D Poisson equation with homogeneous Dirichlet boundary conditions');
    for i=1:size(N,2)
        [e(i), A]=error_PseudoSpectral3D(N(i)-1);
        sizeA(i)=size(A,1);
        sparsity(i)=issparse(A);
        centrosymmetry(i)=iscentrosym(A);
        %condition(i)=cond(full(A));
    end
   Results=table(N',e',sizeA',sparsity',centrosymmetry',...
        'VariableNames',{'number of nodes', 'error',....
        'size of the matrix','sparsity', 'centrosymmetry'})
elseif ecase==4
    N=18:2:22;
    disp('Accuracy of Chebyshev collocation method for 2D diffusion equation with homogeneous Dirichlet boundary conditions');
    for i=1:size(N,2)
        [e(i), A]=error_Poisson_varcoef(N(i)-1,2,10,1);
        sizeA(i)=size(A,1);
        sparsity(i)=issparse(A);
        centrosymmetry(i)=iscentrosym(A);
        condition(i)=cond(full(A));
    end
   Results=table(N',e',sizeA',sparsity',centrosymmetry',condition',...
        'VariableNames',{'number of nodes', 'error',....
        'size of the matrix','sparsity', 'centrosymmetry','condition number'})
elseif ecase==5
    N=8:2:18;
    disp('Accuracy of Chebyshev collocation method for 1D biharmonic equation with homogeneous Dirichlet boundary conditions');
    for i=1:size(N,2)
        [e(i), A]=error_PseudoSpectral1D(N(i),4);
        sizeA(i)=size(A,1);
        sparsity(i)=issparse(A);
        centrosymmetry(i)=iscentrosym(A);
        condition(i)=cond(full(A));
    end
   Results=table(N',e',sizeA',sparsity',centrosymmetry',condition',...
        'VariableNames',{'number of nodes', 'error',....
        'size of the matrix','sparsity', 'centrosymmetry','condition number'})
elseif ecase==6
    N=14:2:26;
    disp('Accuracy of Chebyshev collocation method for 2D biharmonic equation with homogeneous Dirichlet boundary conditions');
    for i=1:size(N,2)
        [e(i), A]=error_Biharmonic_varcoef(N(i)-1,0,1);
        sizeA(i)=size(A,1);
        sparsity(i)=issparse(A);
        centrosymmetry(i)=iscentrosym(A);
        condition(i)=cond(full(A));
    end
    Results=table(N',e',sizeA',sparsity',centrosymmetry',condition',...
        'VariableNames',{'number of nodes', 'error',....
        'size of the matrix','sparsity', 'centrosymmetry','condition number'})
elseif ecase==7
    N=10:2:26;
    disp('Accuracy of Chebyshev collocation method for 2D variable-coefficient biharmonic equation with homogeneous Dirichlet boundary conditions');
    for i=1:size(N,2)
        [e(i), A]=error_Biharmonic_varcoef(N(i)-1,10,1);
        sizeA(i)=size(A,1);
        sparsity(i)=issparse(A);
        centrosymmetry(i)=iscentrosym(A);
        condition(i)=cond(full(A));
    end
    Results=table(N',e',sizeA',sparsity',centrosymmetry',condition',...
        'VariableNames',{'number of nodes', 'error',....
        'size of the matrix','sparsity', 'centrosymmetry','condition number'})
elseif ecase==8
    N=21:1:22;
    disp('Accuracy of Chebyshev collocation method for 3D Helmholtz equation with homogeneous Dirichlet boundary conditions');
    for i=1:size(N,2)
        [e(i), A]=error_Helmholtz3D(N(i)-1,100);
        sizeA(i)=size(A,1);
        sparsity(i)=issparse(A);
        centrosymmetry(i)=iscentrosym(A);
        %condition(i)=cond(full(A));
    end
    Results=table(N',e',sizeA',sparsity',centrosymmetry',...
        'VariableNames',{'number of nodes', 'error',....
        'size of the matrix','sparsity', 'centrosymmetry'})
elseif ecase==9
    N=14:2:26;
    disp('Accuracy of symmetric Legendre collocation method for 2D Poisson equation with homogeneous Dirichlet boundary conditions');
    for i=1:size(N,2)
        [e(i), A]=error_PoissonSymLeg(N(i)-1);
        sizeA(i)=size(A,1);
        sparsity(i)=issparse(A);
        centrosymmetry(i)=iscentrosym(A);
        condition(i)=cond(full(A));
    end
    Results=table(N',e',sizeA',sparsity',centrosymmetry',condition',...
        'VariableNames',{'number of nodes', 'error',....
        'size of the matrix','sparsity', 'centrosymmetry','condition number'})
elseif ecase==10
    N=16:2:20;
    disp('Accuracy of symmetric Legendre collocation method for the PDE \Delta u+u=f(x,y) with Neumann boundary conditions');
    for i=1:size(N,2)
        [e(i), A]=error_Neumann2D(N(i)-1);
        sizeA(i)=size(A,1);
        sparsity(i)=issparse(A);
        centrosymmetry(i)=iscentrosym(A);
        condition(i)=cond(full(A));
    end
    Results=table(N',e',sizeA',sparsity',centrosymmetry',condition',...
        'VariableNames',{'number of nodes', 'error',....
        'size of the matrix','sparsity', 'centrosymmetry','condition number'})
elseif ecase==11
    N=10:2:18;
    disp('Accuracy of symmetric Legendre collocation method for 2D Poisson equation with Robin boundary conditions');
    for i=1:size(N,2)
        [e(i), A]=error_PoissonRobin2D(N(i)-1);
        sizeA(i)=size(A,1);
        sparsity(i)=issparse(A);
        centrosymmetry(i)=iscentrosym(A);
        condition(i)=cond(full(A));
    end
    Results=table(N',e',sizeA',sparsity',centrosymmetry',condition',...
        'VariableNames',{'number of nodes', 'error',....
        'size of the matrix','sparsity', 'centrosymmetry','condition number'})
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
end
xlabel('N')
ylabel('error')
plotformat(1.5,6)

function c=iscentrosym(A)
J=flipud(eye(size(A,1)));B=0.5*(A+J*A*J);
diffnorm=norm(A-B);
if diffnorm<=1e-8
    c=true;
else 
    c=false;
end
end
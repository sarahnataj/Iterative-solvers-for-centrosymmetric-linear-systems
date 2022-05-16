function [D, node]=PSDirv(n,method)
%Calculate the first order spectral differentiation matrix and Guass Lobatto nodes
%method 1=Legendre  2=Chebyshev
n1=n+1; D=zeros(n1);
if method==1 % Legendre
    j=1:n-2; d=.5*sqrt(j.*(j+2)./(j+.5)./(j+1.5));
    node=eig(diag(d,1)+diag(d,-1)); % Legendre Lobatto pts
    node=[-1; flipud(node); 1];
    legn=legendre(n,node); legn=legn(1,:); % L_n(nodes)
    for i=1:n1 % calculate deriv matrix
        legni=legn(i); nodei=node(i);
        for j=1:n1
            if i ~= j
                D(i,j)=legni/legn(j)/(nodei-node(j));
            end
        end
    end
    D(1,1)=n*n1/4; D(n1,n1)=-n*n1/4; w=2/n/(n+1)./legn.^2; w=w';
else % Cheb
    j=0:n; d=ones(n1,1); d(1)=2; d(n1)=2;
    node=cos(pi*j/n)';
    for i=1:n1 % calculate deriv matrix
        nodei=node(i); di=d(i);  D(1,1)=(2*n*n+1)/6; D(n1,n1)=-D(1,1);
        if i>1 && i<n1
            D(i,i)=-nodei/2/(1-nodei*nodei);
        end
        for j=1:n1
            if i ~= j
                D(i,j)=(-1)^(i+j)*di/d(j)/(nodei-node(j));
            end
        end
    end
end
end
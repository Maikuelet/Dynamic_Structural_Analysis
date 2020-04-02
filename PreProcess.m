%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%       PREPROCESS         %%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [dim, T2 ] = PreProcess(input)

%% Dimensions
dim.Ndim = size(input.position,1);      % Number of dimensions (DOFs for each node)
dim.Nnodes = size(input.position,2);    % Number of nodes
dim.NnodesXelement = size(input.T,1);   % Number of nodes for each element
dim.Nelements = size(input.T,2);        % Number of elements
dim.Ndofs = dim.Nnodes*dim.Ndim;        % Total number of degrees of freedom



%% Conectivity table 
%  Defined in terms of the degrees of freedom

nrows=dim.NnodesXelement*dim.Ndim;
T2=zeros(dim.Nelements,nrows);

for i=1:dim.Nelements
    for j=1:nrows
        if (j==1)
            T2(j,i)=input.T(1,i)*2-1;
        end
        if(j==2)
             T2(j,i)=input.T(1,i)*2;
        end
        if(j==3)
           T2(j,i)=input.T(2,i)*2-1;
        end
        if(j==4)
           T2(j,i)=input.T(2,i)*2;  
        end
    end
end

end
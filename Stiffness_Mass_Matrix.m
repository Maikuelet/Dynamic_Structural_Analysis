%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%  STIFFNESS MATRIX  %%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%Function in charge of generating the Global Stiffness Matrix


function [KG,MG]  = Stiffness_Mass_Matrix(dim,input,T2)

% Stiffness and mass matrix for each bar element
Kel = zeros(dim.NnodesXelement*dim.Ndim,dim.NnodesXelement*dim.Ndim,dim.Nelements);
KG = zeros(dim.Ndofs,dim.Ndofs);
MG.C = zeros(dim.Ndofs,dim.Ndofs);
MG.L = zeros(dim.Ndofs,dim.Ndofs);
MG.O = zeros(dim.Ndofs,dim.Ndofs);

for e=1:dim.Nelements
    
    % Element Position
    x1e=input.position(1,input.T(1,e));
    x2e=input.position(1,input.T(2,e));
    y1e=input.position(2,input.T(1,e));
    y2e=input.position(2,input.T(2,e));
    
    % Element Lenght and orientation
    le=sqrt((x2e-x1e)^2 +(y2e-y1e)^2);
    se=((y2e-y1e)/le);                  % sine   of Element
    ce=((x2e-x1e)/le);                  % cosine of Element
    
    
    %Element Matrix 
    %       Ae * Ee   /   le
    Ee   = input.mat(1,input.Mind(1,e));
    Ae   = input.mat(2,input.Mind(1,e));
    rhoe = input.mat(3,input.Mind(1,e));
    ke=((Ae*Ee)/le)...
        *([ce^2 ce*se -ce^2 -ce*se;
           ce*se se^2 -ce*se -se^2;
           -ce^2 -ce*se ce^2 ce*se;
           -ce*se -se^2 ce*se se^2;
        ]);
    
    % Store Element Matrix
    for r=1:dim.Ndim*dim.NnodesXelement
        for s=1:dim.Ndim*dim.NnodesXelement
            Kel(r,s,e)=ke(r,s);
        end   
    end
    
    %% Mass Matrix
    Rot=[ce -se 0 0;
              se  ce 0 0;
              0   0 ce -se;
              0   0 se ce];
    % Consistent
    Mc = ((rhoe*le*Ae)/2)*[2/3   0  1/3   0;
                           0   2/3  0   1/3;
                           1/3  0  2/3   0;
                           0   1/3  0   2/3]; 
    % Lumped                  
    Ml = ((rhoe*le*Ae)/2)*[1 0 0 0;
                           0 1 0 0;
                           0 0 1 0;
                           0 0 0 1]; 
    % Optimal
    Mo = ((rhoe*le*Ae/2))*[5/6  0  1/6  0;
                           0    1  0    0;
                           1/6  0  5/6  0;
                           0    0  0    1]; % optimal    
    Mc = Rot*Mc*Rot';
    Ml = Rot*Ml*Rot';
    Mo = Rot*Mo*Rot';
        
    for r=1:dim.Ndim*dim.NnodesXelement
        for s=1:dim.Ndim*dim.NnodesXelement
            Mcc (r,s,e) = Mc(r,s);
            Mll (r,s,e) = Ml(r,s);
            Moo (r,s,e) = Mo(r,s);
        end   
    end
   
end

% if (a == 1)
%     M_aux = Mcc; % The global matrix will be based on the consistent one
% elseif(a == 2 )
%     M_aux = Mll; % The global matrix will be based on the lumped one
% elseif(a == 3 )
%     M_aux = Moo; % The global matrix will be based on the optimal one
% end


%%  Global Matrix Storage
for e=1:dim.Nelements                             % For each Beam
    for i=1:dim.Ndim*dim.NnodesXelement           % Local degree of freedom
        I=T2(i,e);                                % Global dof
        for j= 1:dim.Ndim*dim.NnodesXelement      % Local dof
            J=T2(j,e);                            % Global dog
            KG(I,J)=KG(I,J)+ Kel(i,j,e);
            MG.C(I,J)= MG.C(I,J)+ Mcc(i,j,e);
            MG.L(I,J)= MG.L(I,J)+ Mll(i,j,e);
            MG.O(I,J)= MG.O(I,J)+ Moo(i,j,e);
        end
    end
end

end
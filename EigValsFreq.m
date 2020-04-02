%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%  EIGENVALUES AND EIGENFREQ  %%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [EigVal,Omega,direct,modal] = EigValsFreq(input,dim,KG,MG,eign)


%% Restricted degrees of freedom
% Vector conaining all the fixed degrees of freedom
NfixNod = size(input.fixNod,2);              
vR = zeros(1,NfixNod);     
fixnodeposi = input.fixNod(1,:)*2-2 + input.fixNod(2,:);   % Node*2dim -2 + DOF
vR(1,fixnodeposi) = input.fixNod(1,:)*2-2 + input.fixNod(2,:);

% Vector conaining all the free degrees of freedom
vL = setdiff(1:dim.Ndofs,vR);

%Reducing the 0's of an array
vL = vL(vL~=0);

%% Reduced system matrices
KLL = KG(vL,vL);
MLL.L = MG.L(vL,vL);
MLL.O = MG.O(vL,vL);
MLL.C = MG.C(vL,vL);

%Extracting the eigenmodes & eigenfrequencies
[EigVal.L,Omega.L] = eig(KLL,MLL.L);
[EigVal.O,Omega.O] = eig(KLL,MLL.O);
[EigVal.C,Omega.C] = eig(KLL,MLL.C);

for i = 1:eign.vals
    min_val = min(abs(EigVal.L(:,i)));
    EigVal.L(i,:) = EigVal.L(i,:)./min_val;
    
    min_val = min(abs(EigVal.O(:,i)));
    EigVal.O(i,:) = EigVal.O(i,:)./min_val;
    
    min_val = min(abs(EigVal.C(:,i)));
    EigVal.C(i,:) = EigVal.C(i,:)./min_val;
end

%% DIRECT METHOD COMPUTATION

for i = 1:length(eign.freq)
    %Receptance matrice assembly
    QL(:,:,i) = -eign.freq(i)^2*MLL.L + KLL;
    QO(:,:,i) = -eign.freq(i)^2*MLL.O + KLL;
    QC(:,:,i) = -eign.freq(i)^2*MLL.C + KLL;
    direct.qL(i,:) = diag(inv(QL(:,:,i)));
    direct.qO(i,:) = diag(inv(QO(:,:,i)));
    direct.qC(i,:) = diag(inv(QC(:,:,i)));
end

%% MODAL MODE COMPUTATION

 modal.qL = zeros(length(eign.freq),size(Omega.L,1));
 modal.qO = zeros(length(eign.freq),size(Omega.O,1));
 modal.qC = zeros(length(eign.freq),size(Omega.C,1));

for i = 1:length(eign.freq)
    modal.qL(i,:) = diag(inv(-eign.freq(i)^2*eye(size(Omega.L,1))+Omega.L));
    modal.qO(i,:) = diag(inv(-eign.freq(i)^2*eye(size(Omega.O,1))+Omega.O));
    modal.qC(i,:) = diag(inv(-eign.freq(i)^2*eye(size(Omega.C,1))+Omega.C));
end

end 
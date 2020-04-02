%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%  STATIC STRUCTURE 2D SOLVER   %%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%     STRUCTURAL DYNAMICS       %%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%-- AUTHORS --%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%  Miquel Altadill Llasat  %%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%     Julia Blanch Sánchez        %%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot and save deformed and undeformed shape of 2D structure

clc
clear all
close all

%% =========================== INPUT DATA ============================== %%
InputData

%% ==========================  PREPROCESS =============================  %%
[dim, T2 ] = PreProcess(input);

%% ============================= SOLVER ===============================  %%

[KG,MG]  = Stiffness_Mass_Matrix(dim,input,T2);

%% ========================= FREQ. RESPONSE ===========================  %%
[EigVal,Omega,direct,modal] = EigValsFreq(input,dim,KG,MG,eign);

%% ========================== POSTPROCESS =============================  %%
postproces(modal,direct,Omega,eign,m1,KG,EigVal);


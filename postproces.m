%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%       POSTPROCESS        %%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function postproces(modal,direct,Omega,eign,m1,KG,EigVal)

%% PLOT 1
h1 = figure(1);
maxfreq = eign.freq(length(eign.freq))/(2*pi);  % In rad/s


for i=1:size(Omega.O,1)
    
    x = 'Frequency (Hz)';
    y = '|u| / F';
    
    subplot (1,3,1)
    ylim([0 0.8e-3])
    xlim([0 maxfreq])
    plot(eign.freq./(2*pi),(abs(direct.qL(:,i))));
    title('Lumped','Interpreter','latex');
    xlabel(x ,'Interpreter','latex');
    ylabel(y,'Interpreter','latex');
    hold on
    
    subplot (1,3,2)
    ylim([0 0.8e-3])
    xlim([0 maxfreq])
    plot(eign.freq./(2*pi),(abs(direct.qO(:,i))));
    title('Optimal','Interpreter','latex');
    xlabel(x ,'Interpreter','latex');
    ylabel(y,'Interpreter','latex');
    hold on

    subplot (1,3,3)
    ylim([0 0.8e-3])
    xlim([0 maxfreq])
    plot(eign.freq./(2*pi),(abs(direct.qC(:,i))));
    legend('Receptance 1','Receptance 2');
    title('Consistent','Interpreter','latex');
    xlabel(x ,'Interpreter','latex');
    ylabel(y,'Interpreter','latex');
    hold on
end

%% PLOT 2
h2 = figure(2);
for i=1:size(Omega.O,1)
    x = 'Frequency (Hz)';
    y = {'$Q / \widetilde{F}$'};
    
    subplot (1,3,1)
    ylim([0 1e-4])
    xlim([0 maxfreq])
    title('Lumped','Interpreter','latex');
    plot(eign.freq./(2*pi),(abs(modal.qL(:,i))));
    xlabel(x,'Interpreter','latex');
    ylabel(y,'Interpreter','latex');
    hold on
    
    subplot (1,3,2)
    ylim([0 1e-4])
    xlim([0 maxfreq])
    plot(eign.freq./(2*pi),(abs(modal.qO(:,i))));
    title('Optimal','Interpreter','latex');
    xlabel(x,'Interpreter','latex');
    ylabel(y,'Interpreter','latex');
    hold on
    
    subplot (1,3,3)
    ylim([0 1e-4])
    xlim([0 maxfreq])
    plot(eign.freq./(2*pi),(abs(modal.qC(:,i))));
    title('Consistent','Interpreter','latex');
    xlabel(x,'Interpreter','latex');
    ylabel(y,'Interpreter','latex');
    legend('Receptance 1','Receptance 2');
    hold on
end


%% Results for displaying
Omega.O = sqrt((diag(Omega.O)));
freqO=Omega.O/(2*pi);
Omega.C = sqrt((diag(Omega.C)));
freqC=Omega.C/(2*pi);
Omega.L= sqrt((diag(Omega.L)));
freqL=Omega.L/(2*pi);

%% Displays 
disp('-----------------------INFORMATION-------------------------');
disp('EigFreq Optimal [Hz]:');
disp(freqO);
disp('EigFreq Lumped [Hz]:');
disp(freqL);
disp('EigFreq Consistent [Hz]:');
disp(freqC);
disp('-------------------------------------------------------------');
disp('EigVal Optimal [Hz]:');
disp(EigVal.O);
disp('EigVal Lumped [Hz]:');
disp(EigVal.L);
disp('EigVal Consistent [Hz]:');
disp(EigVal.C);
disp('-------------------------------------------------------------');
disp('Stiffness Matrix [Hz]:');
disp(KG);


%% Plot Storage
pathh     = pwd;
myfolder = 'Plots';
f1 = fullfile(pathh , myfolder);
mkdir(f1);       % Creates Folder

f = fullfile(f1 , sprintf('Direct_Rho%d.png', m1.rho));
saveas(h1,f);
f = fullfile(f1 , sprintf('Modal_Rho%d.png', m1.rho));
saveas(h2,f);


end
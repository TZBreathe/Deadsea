%% Check Var(dQ) correlation within present datasets up to only 100 cycle
clear;
addpath(genpath('G:\09_Research_Development_Projects\00_Degradation data sets'));

%Temp lazy implementation

VardQ=struct('Elbrus',[],'K2',[],'Varta',[],'VartaBr',[]);
Cap100=struct('Elbrus',[],'K2',[],'Varta',[],'VartaBr',[]);


feature_extraction_Elbrus;
VardQ.Elbrus=VarQ;
Cap100.Elbrus=EndCap;

clear EndCap VarQ;

feature_extraction_K2;
VardQ.K2=VarQ;
Cap100.K2=EndCap;

clear EndCap VarQ;

feature_extraction_varta;
VardQ.Varta=VarQ;
Cap100.Varta=EndCap;

clear EndCap VarQ;

feature_extraction_vartaBr;
VardQ.VartaBr=VarQ;
Cap100.VartaBr=EndCap;

% Var_dQ=[VardQ{1} VardQ{2} VardQ{3} VardQ{4}];
% Cap_100=[Cap100{1} Cap100{2} Cap100{3} Cap100{4}];

plot(VardQ.Elbrus,Cap100.Elbrus, 'o'); hold on;
plot(VardQ.K2,Cap100.K2, 's');
plot(VardQ.Varta,Cap100.Varta, 'x');
plot(VardQ.VartaBr,Cap100.VartaBr, '*');
legend('Elbrus cycle 150','K2 cycle 600','Varta client cycle 350','Varta Br cycle 100','Location','southeast');
hold off;
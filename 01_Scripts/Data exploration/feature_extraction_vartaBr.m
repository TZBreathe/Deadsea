%% starter code for extracting features for degradation regression model

% Varta Br data more than 100 cycles
addpath(genpath('G:\09_Research_Development_Projects\00_Degradation data sets'));
%
runNumbers=[740 1018]; 
Br{1}=appendDegData(runNumbers);
runNumbers=[741 1019]; 
Br{2}=appendDegData(runNumbers);
runNumbers=[742 1020]; 
Br{3}=appendDegData(runNumbers);
runNumbers=[732 1013]; 
Br{4}=appendDegData(runNumbers);
runNumbers=[733 1014]; 
Br{5}=appendDegData(runNumbers);
runNumbers=[728 1012]; 
Br{6}=appendDegData(runNumbers);
runNumbers=[735 1016]; 
Br{7}=appendDegData(runNumbers);
runNumbers=[736 1017]; 
Br{8}=appendDegData(runNumbers);
Br{9}=gdFun.Load_Breathe_Run(901);
Br{10}=gdFun.Load_Breathe_Run(902);
Br{11}=gdFun.Load_Breathe_Run(903);
Br{12}=gdFun.Load_Breathe_Run(898);
Br{13}=gdFun.Load_Breathe_Run(899);
Br{14}=gdFun.Load_Breathe_Run(900);
Br{15}=gdFun.Load_Breathe_Run(895);
Br{16}=gdFun.Load_Breathe_Run(896);
Br{17}=gdFun.Load_Breathe_Run(897);

MaxCap=3;

%% 
%variance of dQ(V)

for i=1:length(Br)
    
    
    ind_10=find([Br{i}.RunData.dataTable{:,'cycleNumber'}==5]&(Br{i}.RunData.dataTable{:,'currCell'}<0));
    ind_50=find([Br{i}.RunData.dataTable{:,'cycleNumber'}==51]&(Br{i}.RunData.dataTable{:,'currCell'}<0));
    
    VQ_10_tmp=[Br{i}.RunData.dataTable{ind_10,'voltCell'} Br{i}.RunData.dataTable{ind_10,'ahDchrge'}];
      [~,ia,ic]=unique(VQ_10_tmp(:,1));   %remove repeated V values
    VQ_10=VQ_10_tmp(ia,:);
      VQ_10(:,2)=VQ_10(:,2)./MaxCap; %Normalise capacity, make sense or not?
    
    V_interp=linspace(VQ_10(1,1),VQ_10(end,1),1000);
    Q_interp=interp1(VQ_10(:,1),VQ_10(:,2),V_interp);
    VQ_10=[V_interp; Q_interp];
    
    VQ_50_tmp=[Br{i}.RunData.dataTable{ind_50,'voltCell'} Br{i}.RunData.dataTable{ind_50,'ahDchrge'}];
     [~,ia,ic]=unique(VQ_50_tmp(:,1));   %remove repeated V values
     VQ_50=VQ_50_tmp(ia,:);
    VQ_50(:,2)=VQ_50(:,2)./MaxCap;
    
    V_interp=linspace(VQ_50(1,1),VQ_50(end,1),1000);
    Q_interp=interp1(VQ_50(:,1),VQ_50(:,2),V_interp);
    VQ_50=[V_interp; Q_interp];
    
     dQ=abs(VQ_10(2,:)-VQ_50(2,:));
     VarQ(i)=var(dQ);
     AvgQ(i)=mean(dQ);
     MaxQ(i)=max(dQ);
     
       EndCap(i)=Br{i}.RunData.cycleTable{100,'ahDchrge'}./MaxCap; 
     %     plot(VQ_10(1,:),VQ_10(2,:),VQ_50(1,:),VQ_50(2,:));
%        plot(dQ,VQ_4(1,:));

end
     
% 
% plot(EndCap, VarQ,'O');
% figure
% plot(EndCap, AvgQ,'O');
% figure
% plot(EndCap, MaxQ,'O');

%%
function [dataoutput]=appendDegData(runNumbers)

    br=gdFun.Load_Breathe_Run(runNumbers); 

    	for i_r = 1:length(br)
		
        %%removing the tcCounter3 column
        if any(contains(br{i_r}.RunData.dataTable.Properties.VariableNames,'tcCounter3'))
            col2del = find(contains(br{i_r}.RunData.dataTable.Properties.VariableNames,'tcCounter3'));
            br{i_r}.RunData.dataTable(:,col2del) = [];
        end
            
		if i_r ==1 
			dtCell = br{i_r}.RunData.dataTable;
			ctCell = br{i_r}.RunData.cycleTable;
		else
			dtCell = vertcat(dtCell,br{i_r}.RunData.dataTable);
			
			% Extend the table to add the keys from the new table (not all tables have the same keys)
			tmpTableAdd = br{i_r}.RunData.cycleTable;
			tmpTableAdd(:,:) = []; % Make an empty table with only headers
			ctExtended = outerjoin(tmpTableAdd,ctCell,'MergeKeys',true);
			
			ctCell = vertcat(ctExtended,br{i_r}.RunData.cycleTable);
		end
	end

	dataoutput.RunData.cycleTable = ctCell;
	dataoutput.RunData.dataTable = dtCell;

end
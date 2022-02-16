% Extracts discharge curves in a degradation dataset and re-sample them to same length.
% Also extracts capacities

function [Vcell, yS] = get_disch_curves(runNo)

    data=gdFun.Load_Breathe_Run(runNo);
    %load capacities
    i=0;
    y=[];
    for j=1:length(data)
        for i=1:height(data{j}.RunData.cycleTable)     
            y=vertcat(y,data{j}.RunData.cycleTable{i,'ahDchrge'});
        end
    end 
    % yS=smooth(y,5)'; %smoothing
    yS=y';
    
    % Extract discharge curves and resample them
    ii=0;
    Vtmp={};
    Vcell={};
    for j=1:length(data)
   
        Idx_disch=find((data{1,j}.RunData.dataTable{:,'currCell'}<0));
        V_tmp=data{1,j}.RunData.dataTable{Idx_disch,'voltCell'};
        Idx_bp=find(diff(V_tmp)>0.8);
        %first disch curve
        Vtmp{1}=V_tmp(1:Idx_bp(1));
        ii=ii+1;
         Capcell{ii}=(1:length(Vtmp{1}))./length(Vtmp{1});Vcell{ii}=Vtmp{1};
        %second to second-to-last
        for i=2:length(Idx_bp)
            Vtmp{i}=V_tmp(Idx_bp(i-1)+1:Idx_bp(i));
            Captmp{i}=(1:length(Vtmp{i}))./length(Vtmp{i});
            ii=ii+1;
            Vcell{ii}= Vtmp{i};
            Capcell{ii}= Captmp{i};
    %         if length(Vcell{ii+j})>40000 % remove the slow disch cycles
    %             Vcell{ii+j}=[];
    %             Capcell{ii+j}=[];
    %         end      
        end
        %last discharge curve
         Vtmp{i+1}=V_tmp(Idx_bp(i)+1:end);
         Captmp{i+1}=(1:length(Vtmp{i+1}))./length(Vtmp{i+1});
         ii=ii+1;
         Vcell{ii}= Vtmp{i+1};
         Capcell{ii}= Captmp{i+1};   
    end
    
    for i=1:length(Vcell)
         if length(Vcell{i})<1000
             Vcell{i}=[];
             Capcell{i}=[];
         end
    end
    % shouldn't need the following lines if data loaded correctly
    Vcell =Vcell(~cellfun('isempty',Vcell)) ; 
    Capcell =Capcell(~cellfun('isempty',Capcell)) ; 
    % resample so each disch curve has the same length (500)
    for i=1:length(Vcell)
        Vcell_intp{i}=interp1(Capcell{i},Vcell{i},linspace(Capcell{i}(1),1,500));
    end
     
     Vcell=[Vcell_intp{:}];   
end
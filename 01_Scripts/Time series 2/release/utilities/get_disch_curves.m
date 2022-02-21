% Extracts discharge curves in a degradation dataset and re-sample them to same length.
% Also extracts capacities

function [Vcell, yS, IR] = get_disch_curves(runNo)

    data=gdFun.Load_Multiple_Runs(runNo,true);
    %load capacities
    y(1)=data.RunData.cycleTable{1,'ahDchrge'};
    for i=2:height(data.RunData.cycleTable)    
        y(i)=data.RunData.cycleTable{i,'ahDchrge'};
        if y(i)<y(i-1)*0.8    % replace oddities in data where capacity suddenly drops for certain cycles
           y(i)= [];
        end
    end
    cell_cap=y(2);
    y=y./cell_cap;   
    yS=smooth(y,'rlowess');
    
    % Extract discharge curves and resample them
    iVtmp={};
    Capcell={};
    Vcell={};
    Idx_disch=find((data.RunData.dataTable{:,'currCell'}<0));
    V_tmp=data.RunData.dataTable{Idx_disch,'voltCell'};
    Idx_bp=find(diff(V_tmp)>0.7);
    %first disch curve
    Vtmp{1}=V_tmp(1:Idx_bp(1)); Capcell{1}=(1:length(Vtmp{1}))./length(Vtmp{1});Vcell{1}=Vtmp{1};
    %second to second-to-last
    for i=2:length(Idx_bp)
        Vtmp{i}=V_tmp(Idx_bp(i-1)+1:Idx_bp(i));
        Captmp{i}=(1:length(Vtmp{i}))./length(Vtmp{i});
        Vcell{i}= Vtmp{i};
        Capcell{i}= Captmp{i};
    end
    %last discharge curve
    Vtmp{i+1}=V_tmp(Idx_bp(i)+1:end);
    Captmp{i+1}=(1:length(Vtmp{i+1}))./length(Vtmp{i+1}); 
    Vcell{i+1}= Vtmp{i+1};
    Capcell{i+1}= Captmp{i+1};   
    %remove partial/empty cycles
    for i=1:length(Vcell)
         if length(Vcell{i})<1000
             Vcell{i}=[];
             Capcell{i}=[];
         end
    end
    Vcell =Vcell(~cellfun('isempty',Vcell)) ; 
    Capcell =Capcell(~cellfun('isempty',Capcell)) ; 
    % resample so each disch curve has the same length (500)
    for i=1:length(Vcell)
        Vcell_intp{i}=interp1(Capcell{i},Vcell{i},linspace(Capcell{i}(1),1,500));
    end
     
     Vcell=[Vcell_intp{:}];   
     
     %get IR
     idx=[];
     idx(1)=1;
     for i=2:length(Idx_disch)
        if Idx_disch(i)-Idx_disch(i-1)>100
            idx=vertcat(idx,i);
        end
    end
    for i=1:length(idx)
        IR(i)=100*data.RunData.dataTable{Idx_disch(idx(i))-1,'voltCell'}-data.RunData.dataTable{Idx_disch(idx(i)),'voltCell'};
    end
    IR=rmoutliers(IR,'movmean',5);
    if IR(1)-IR(2)>3*(IR(2)-IR(3))
        IR(1)=[];
    end
end
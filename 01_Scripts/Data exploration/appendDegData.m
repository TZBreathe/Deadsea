function [dataoutput]=appendDegData(runNumbers)

    br=gdFun.Load_Breathe_Run(runNumbers); 

    	for i_r = 1:length(br)
		
        %%removing the tcCounter3 column
        if any(contains(br{i_r}.RunData.dataTable.Properties.VariableNames,'tcCounter3'))
            col2del = find(contains(br{i_r}.RunData.dataTable.Properties.VariableNames,'tcCounter3'));
            br{i_r}.RunData.dataTable(:,col2del) = [];
        end
            
        
        if any(contains(br{i_r}.RunData.dataTable.Properties.VariableNames,'socKf'))
            col2del = find(contains(br{i_r}.RunData.dataTable.Properties.VariableNames,'socKf'));
            br{i_r}.RunData.dataTable(:,col2del) = [];
        end
        
          if any(contains(br{i_r}.RunData.cycleTable.Properties.VariableNames,'timeToCvStart'))
            col2del = find(contains(br{i_r}.RunData.cycleTable.Properties.VariableNames,'timeToCvStart'));
            br{i_r}.RunData.cycleTable(:,col2del) = [];
        end
     
           if any(contains(br{i_r}.RunData.cycleTable.Properties.VariableNames,'timeToCvEnd'))
            col2del = find(contains(br{i_r}.RunData.cycleTable.Properties.VariableNames,'timeToCvEnd'));
            br{i_r}.RunData.cycleTable(:,col2del) = [];
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
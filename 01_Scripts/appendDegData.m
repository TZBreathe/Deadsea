function [dataoutput]=appendDegData(runNumbers)

    br=gdFun.Load_Breathe_Run(runNumbers); 
    if length(br) == 1
        dtCell = br(1).RunData.dataTable;
        ctCell = br(1).RunData.cycleTable;
        startCycle = dtCell.cycleNumber(end);
        startTime = dtCell.timeTest(end);
    else
        for i_r = 1:length(br)
            if i_r ==1 
                dtCell = br{1,i_r}.RunData.dataTable;
                ctCell = br{1,i_r}.RunData.cycleTable;
                startCycle = dtCell.cycleNumber(end);
                startTime = dtCell.timeTest(end);
            else
                br{1,i_r}.RunData.cycleTable.cycleNumber = br{1,i_r}.RunData.cycleTable.cycleNumber + startCycle;
                br{1,i_r}.RunData.dataTable.timeTest = br{1,i_r}.RunData.dataTable.timeTest + startTime;
                br{1,i_r}.RunData.dataTable.cycleNumber = br{1,i_r}.RunData.dataTable.cycleNumber + startCycle;
                
                tmpdataAdd = br{1,i_r}.RunData.dataTable;
                tmpdataAdd(:,:) = []; % Make an empty table with only headers
                dtExtended = outerjoin(tmpdataAdd,dtCell,'MergeKeys',true);
                dtCell = vertcat(dtExtended,br{1,i_r}.RunData.dataTable);
            
			
                % Extend the table to add the keys from the new table (not all tables have the same keys)
                tmpTableAdd = br{1,i_r}.RunData.cycleTable;
                tmpTableAdd(:,:) = []; % Make an empty table with only headers
                ctExtended = outerjoin(tmpTableAdd,ctCell,'MergeKeys',true);
			
                ctCell = vertcat(ctExtended,br{1,i_r}.RunData.cycleTable);
            
                startCycle = dtCell.cycleNumber(end);
                startTime = dtCell.timeTest(end);
            end
        end
    end

	dataoutput.RunData.cycleTable = ctCell;
	dataoutput.RunData.dataTable = dtCell;

end
function [y]=deg_extrap(runNo,pred_len,varargin)

%add path
    folder = fileparts(which(mfilename)); 
    addpath(genpath(folder));
%parse input
    p=inputParser;
    defaultMethod='linear';
    validMethod={'linear','arima','arimax','stl'};
    checkMethod=@(x) any(validatestring(x,validMethod));
    
    defaultPlot=false;
    validPlot={'true','false'};
    checkPlot=@(x) any(validatestring(x,validPlot));
% when train_len==0, the actual extrapolation functions will use whole deg series to train    
    defaultTrain_len=0;
    
    addRequired(p,'runNo',@isnumeric);
    addRequired(p,'pred_len',@isnumeric);
    addParameter(p,'method',defaultMethod,checkMethod);
    addParameter(p,'train_len',defaultTrain_len,@isnumeric);
    addParameter(p,'plot',defaultPlot,checkPlot);
    
    parse(p,runNo,pred_len,varargin{:})
    
    method=p.Results.method;
    train_len=p.Results.train_len;
    plot_test=p.Results.plot;
    
    switch method
        case 'linear'
            %disp('linear');
            y=run_linear(runNo,pred_len,plot_test);
        case 'arima'
            %disp('arima');
            y=run_arima(runNo,pred_len,plot_test,train_len);
        case 'arimax'
            y=run_arimax(runNo,pred_len,plot_test,train_len);
        case 'stl'
            %disp('stl');
            y=run_STL(runNo,pred_len,train_len);
    end
    
end  
    
    
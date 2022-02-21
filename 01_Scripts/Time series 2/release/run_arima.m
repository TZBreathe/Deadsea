function [prediction]=run_arima(run_No,pred_len,varargin)
%add path
    folder = fileparts(which(mfilename)); 
    addpath(genpath(folder));
    
    Br{1}=gdFun.Load_Multiple_Runs(run_No,false);
    n=1;% sampling freq of degradation data

    % get degradation trend cuve
    y(1)=Br{1}.cycleTable{2,'ahDchrge'};
        for j=2:height(Br{1}.cycleTable)/n      
            y(j)=Br{1}.cycleTable{n*j,'ahDchrge'};
             if y(j)<y(j-1)*0.7    % replace oddities in data where capacity suddenly drops for certain cycles
                y(j)= y(j-1);
             end
        end
    y(1)=Br{1}.cycleTable{2,'ahDchrge'};% first cycle capacity data is ususally wrong
    y=y./cell_cap;   
    yS=smooth(y,5,'lowess'); %smoothing
    
    if nargin<4
        train_len=length(y);
    elseif varargin{2}==0
        train_len=length(y);
    else
        train_len=floor(varargin{2}/n);
    end
         
    if nargin<3
        plot_test=false;
    else
        plot_test=varargin{1};
    end
          
    
    %convert to python data type and define model orders to train & compare
    X=py.numpy.array(yS(1:train_len));
    prediction=py.arima_functions.run_arima_auto(X,pred_len);
    prediction=double(prediction);
    
%     if plot_pred==1
%         figure()
%         hold on;
%         plot(prediction(1,:)*n,prediction(2,:),'r--');
%         plot(1:n:n*length(y),y,'bs');
%         legend('forecast','data');
%         xlabel('Cycle');
%         ylabel('Normalised capacity');
%         hold off;
%     end
    
    if plot_test=='true'
        X2=py.numpy.array(yS(1:length(yS)*0.75));
        test=py.arima_functions.run_arima_auto(X,length(yS*0.25)); 
        test=double(test);
        figure()
        hold on;
        plot(test(1,:)*n,test(2,:),'r--','LineWidth',1.5);
        plot(1:n:n*length(y),y,'bs');
        legend('Prediction on test data','data');
        xlabel('Cycle');
        ylabel('Normalised capacity');
        hold off;
    end
end
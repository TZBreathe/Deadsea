function [prediction,best_order]=run_arima(run_No,pred_len,varargin)

    Br{1}=gdFun.Load_Multiple_Runs(run_No,false);
    n=3;% sampling freq of degradation data

    % get degradation trend cuve
    y(1)=Br{1}.cycleTable{2,'ahDchrge'};
        for j=2:height(Br{1}.cycleTable)/n      
            y(j)=Br{1}.cycleTable{n*j,'ahDchrge'};
             if y(j)<y(j-1)/1.5    % replace oddities in data where capacity suddenly drops for certain cycles
                y(j)= y(j-1);
             end
        end
    y(1)=Br{1}.cycleTable{2,'ahDchrge'};% first cycle capacity data is ususally wrong
    
    if nargin<6
       cell_cap=y(1);
    else
        cell_cap=varargin{4};
    end
    
    if nargin<5
        plot_test=0;
    else
        plot_test=varargin{3};
    end
    
    if nargin<4
        plot_pred=0;
    else
        plot_pred=varargin{2};
    end
    
    if nargin<3
        start_cycle=length(y);
    else
        start_cycle=floor(varargin{1}/3);
    end
    
    y=y./cell_cap;   
    yS=smooth(y,5); %smoothing

    %convert to python data type and define model orders to train & compare
    X=py.numpy.array(yS(1:start_cycle));
    cycles=py.numpy.array(1:start_cycle);
    p_values = py.list([1:4]);
    d_values = py.list([1:2]);
    q_values = py.list([1:3]);
    pred_len=int16(pred_len/n);
    
    %find best model order based on rmse of prediction vs test data
    best_order=py.arima_functions.grid_search_orders(X,cycles,p_values,d_values,q_values);
    cycles=py.numpy.array(1:start_cycle+double(pred_len)); %there shouldn't be +1, fix this in py function!
    %use the best model and all capacity data to make future predictions
    prediction=py.arima_functions.run_arima_prediction(X,cycles,pred_len,best_order);
    prediction=double(prediction);
    
    if plot_pred==1
        figure()
        hold on;
        plot(prediction(1,:)*n,prediction(2,:),'r--');
        plot(1:n:n*length(y),y,'bs');
        legend('forecast','data');
        xlabel('Cycle');
        ylabel('Normalised capacity');
        hold off;
    end
    
    if plot_test==1
        test=py.arima_functions.run_arima_model(X,pred_len,best_order); 
        test=double(test);
        figure()
        hold on;
        plot(test(1,:)*n,test(2,:),'r--');
        plot(1:n:n*length(y),y,'bs');
        legend('Prediction on test data','data');
        xlabel('Cycle');
        ylabel('Normalised capacity');
        hold off;
    end
end
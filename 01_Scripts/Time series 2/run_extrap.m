function [y,pred_l,pred_p]=run_extrap(run_No,fit_len,pred_len,varargin)


    n=1;% sampling freq of degradation data
    Br{1}=gdFun.Load_Multiple_Runs(run_No,false);


    % get degradation trend cuve
    y(1)=Br{1}.RunData.cycleTable{1,'ahDchrge'};
        for j=2:height(Br{1}.RunData.cycleTable)/n      
            y(j)=Br{1}.RunData.cycleTable{n*j,'ahDchrge'};
              if y(j)<y(j-1)/1.5    % replace oddities in data where capacity suddenly drops for certain cycles
                 y(j)= y(j-1);
             end
        end
    y(1)=Br{1}.RunData.cycleTable{3,'ahDchrge'}; % First capacity is usually wrong
    
    if nargin<6
       cell_cap=y(1);
    else
        cell_cap=varargin{3};
    end
    
    if nargin<5
        plot_test=0;
    else
        plot_test=varargin{2};
    end
    
    if nargin<4
        plot_pred=0;
    else
        plot_pred=varargin{1};
    end
    
   
    y=y/cell_cap;   
    x=1:length(y);
    
   if fit_len=='all'
        fit_len=length(x)-1;
   end
    
    %linear fit using latest fit_len points
    linear_fit=polyfit(x(end-fit_len:end),y(end-fit_len:end),1);
    y_lfit=polyval(linear_fit,x(end-fit_len:end));
%   err=sum(abs(y_lfit-y(end-fit_len:end)));
    %extrapolate into future pred_len points
    y_lpred=polyval(linear_fit,length(x):(length(x)+pred_len));
    pred_l=[length(x):(length(x)+pred_len);y_lpred];
    
    %power law fit using last fit_len points
    p = polyfit(log(x(end-fit_len:end)),log(y(end-fit_len:end)),1); 
    m = p(1);
    b = exp(p(2));
    y_pfit=b*(x(end-fit_len):x(end)).^m;
     %extrapolate into future pred_len points
    y_ppred=b*(length(x):(length(x)+pred_len)).^m;
    pred_p=[length(x):(length(x)+pred_len);y_ppred];
    
    if plot_pred==1
        figure()
        hold on;
        plot(length(x):(length(x)+pred_len),y_lpred,'r--');
        plot(length(x):(length(x)+pred_len),y_ppred,'g--');
        plot(1:length(y),y,'bs');
        legend('linear extrapoltion','power law extrapolation','data');
        xlabel('Cycle');
        ylabel('Normalised capacity');
        hold off;
    end
    
    if plot_test==1    
        figure()
        hold on;
        plot(x(end-fit_len:end),y_lfit,'r--');
        plot(x(end-fit_len):x(end),y_pfit,'g--');
        plot(1:length(y),y,'bs');
        legend('Linear extrapolation','power law extrapolation','data');
        xlabel('Cycle');
        ylabel('Normalised capacity');
        hold off;
    end
end
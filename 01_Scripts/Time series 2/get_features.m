function [Vmean,Vsa,Vsv]=get_features(Vcell)

    Vmean=[];
    Vsa=[];
    % seasonal decomposition and featurization. Each disch curve is 500
    % data points
    X=py.numpy.array(Vcell);
    X_decomp=double(py.run_tseries.decomp_tseries(X));
%     Vseason=ones(500,length(X_decomp)/500);
    for i=0:length(X_decomp)/500-1
        Vmean(i+1)=mean(X_decomp(1,500*i+1:500*(i+1)));    
        Vseason(:,i+1)=X_decomp(2,500*i+1:500*(i+1));
        Vsa(i+1)=trapz(abs(Vseason(:,i+1)-Vseason(:,1)));
        Vsv(i+1)=var(abs(Vseason(:,i+1)-Vseason(:,1)))*100;
    end
end
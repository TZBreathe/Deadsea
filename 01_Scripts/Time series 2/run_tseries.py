# -*- coding: utf-8 -*-
"""
Created on Wed Dec  1 11:20:24 2021

@author: Teng
"""

import numpy as np
import pandas as pd
from matplotlib import pyplot
from statsmodels.tsa.seasonal import seasonal_decompose
from statsmodels.tsa.seasonal import STL
from numpy import loadtxt
from statsmodels.tsa.arima.model import ARIMA
from sklearn.metrics import mean_squared_error
import warnings
warnings.filterwarnings("ignore")

#seasonal decomposition
def decomp_tseries(X):
    #result=seasonal_decompose(X,model='additive',period=500,extrapolate_trend='freq')
    stl=STL(X,period=500,robust=True)
    result=stl.fit()
    return np.vstack((result.trend,result.seasonal))


#Evaluate rmse of a given arima model order
def evaluate_arima_model(X,arima_order):
    
    train_size = int(len(X) * 0.75)
    train, test = X[0:train_size], X[train_size:]
    
    #fit and test
    model = ARIMA(train, order=(arima_order))
    model_fit = model.fit()
    yhat = model_fit.predict(len(train),len(X)-1)
    
    # calculate out of sample error
    error = mean_squared_error(test, yhat)
    return error

#grid search model model orders defined below
def grid_search_orders(dataset,p_values,d_values,q_values):
    import numpy as np
    import pandas as pd
    from matplotlib import pyplot
    from numpy import loadtxt
    from statsmodels.tsa.arima.model import ARIMA
    from sklearn.metrics import mean_squared_error
    import warnings
    warnings.filterwarnings("ignore")
    dataset = dataset.astype('float32')
    best_score, best_cfg = float("inf"), None
    for p in p_values:
        for d in d_values:
           for q in q_values:
              order = (p,d,q)
              try:
                  rmse = 100*evaluate_arima_model(dataset, order)
                  if rmse < best_score:
                      best_score, best_cfg = rmse, order
                      #print('ARIMA%s RMSE=%.3f' % (order,rmse))                    
              except:
                  continue
    #print('Best ARIMA%s RMSE=%.3f' % (best_cfg, best_score))
    return best_cfg

#Prediction using arima
def run_arima_prediction(X,prediction_len,arima_order):
        
    train_start=0
    train= X[train_start:]
    pred=np.zeros(prediction_len)
    
    #run prediction  
    model = ARIMA(train, order=(arima_order))
    model_fit = model.fit()
         
    pred= model_fit.predict(len(train),len(train)+prediction_len-1)
    X_pred=np.arange(len(X),len(X)+prediction_len)
    return np.vstack((X_pred,pred))

#run arima, same as evaluate arima
def run_arima_model(X,prediction_len,arima_order):
    train_size = int(len(X) * 0.75)
    train, test = X[0:train_size], X[train_size:]
    
    #fit and test
    model = ARIMA(train, order=(arima_order))
    model_fit = model.fit()
    yhat = model_fit.predict(len(train),len(X)-1)
    xhat=np.arange(len(train),len(X))
    return np.vstack((xhat,yhat))

# -*- coding: utf-8 -*-
"""
Created on Wed Dec  1 11:20:24 2021

@author: Teng
"""

import numpy as np
import pandas as pd
from matplotlib import pyplot
from numpy import loadtxt
from statsmodels.tsa.arima.model import ARIMA
from statsmodels.tsa.statespace.sarimax import SARIMAX
from sklearn.metrics import mean_squared_error
from statsmodels.tsa.seasonal import seasonal_decompose
from statsmodels.tsa.seasonal import STL
import warnings
warnings.filterwarnings("ignore")

#seasonal decomposition
def decomp_tseries(X):
    #result=seasonal_decompose(X,model='additive',period=500,extrapolate_trend='freq')
    stl=STL(X,period=500,robust=True)
    result=stl.fit()
    return np.vstack((result.trend,result.seasonal))

#Evaluate rmse of a given arima model order
def evaluate_arima_model(X,x2,arima_order):
    
    train_size = int(len(X) * 0.75)
    train, test = X[0:train_size], X[train_size:]
    ex=x2[0:train_size]
    #fit and test
    model=SARIMAX(train,exog=ex,order=(arima_order))
    model_fit = model.fit()
    ex_p=x2[train_size:len(x2)]
    yhat=model_fit.predict(len(train),len(x2)-1,exog=[ex_p])
    # calculate out of sample error
    error = mean_squared_error(test, yhat)
    return error

#grid search model model orders defined below
def grid_search_orders(dataset,x2,p_values,d_values,q_values):
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
                  rmse = 100*evaluate_arima_model(dataset,x2,order)
                  if rmse < best_score:
                      best_score, best_cfg = rmse, order
                      #print('ARIMA%s RMSE=%.3f' % (order,rmse))                    
              except:
                  continue
    #print('Best ARIMA%s RMSE=%.3f' % (best_cfg, best_score))
    return best_cfg

#Prediction using arima
def run_arima_prediction(X,x2,prediction_len,arima_order):
    #x2 vector of cycle number up to extrapolated cycle
    train= X
    ex=x2[0:len(X)]
    pred=np.zeros(prediction_len)
    ex_p=x2[len(X):len(x2)]
    
    #run prediction  
    model=SARIMAX(train,exog=ex,order=(arima_order))
    model_fit = model.fit()   
    pred= model_fit.predict(len(X),len(x2)-1,exog=[ex_p])
    X_pred=np.arange(len(X),len(X)+prediction_len)
    return np.vstack((X_pred,pred))

'''
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
'''

# run arima with ex factors
def run_arima(x1,x2,prediction_len,arima_order):
    from statsmodels.tsa.statespace.sarimax import SARIMAX
    train_size = int(len(x1) * 0.75)
    train, test = x1[0:train_size], x1[train_size:]
    ex=x2[0:train_size]
    model=SARIMAX(train,exog=ex,order=(arima_order))
    model_fit = model.fit()
    ex_p=x2[train_size:len(x2)]
    yhat=model_fit.predict(len(train),len(x2)-1,exog=[ex_p])
    xhat=np.arange(len(train),len(x1))
    return np.vstack((xhat,yhat))

#autoarima test
def run_arima_auto(x1,prediction_len):
    from pmdarima.arima import auto_arima
    train_size = int(len(x1) * 0.75)
    train, test = x1[0:train_size], x1[train_size:]
    model_fit=auto_arima(train,error_action='ignore',supress_warnings=True,maxiter=5)
    yhat=model_fit.predict(n_periods=prediction_len)
    xhat=np.arange(len(train),len(train)+prediction_len)
    return np.vstack((xhat,yhat))  

#autoarima test with ex
def run_arima_auto_ex(x1,x2,prediction_len):
    from pmdarima.arima import auto_arima
    train_size = int(len(x1) * 0.75)
    train, test = x1[0:train_size], x1[train_size:]
    ex=x2[0:train_size]
    model_fit=auto_arima(train,ex,error_action='ignore',supress_warnings=True)
    ex_p=x2[train_size:len(x2)]
    yhat=model_fit.predict(n_periods=prediction_len,ex_p)
    xhat=np.arange(len(train),len(train)+prediction_len)
    return np.vstack((xhat,yhat))  

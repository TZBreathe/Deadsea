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
from statsmodels.tsa.holtwinters import ExponentialSmoothing as HWES
from statsmodels.tsa.seasonal import STL
import warnings
warnings.filterwarnings("ignore")


#HW smoothing trial
def run_hw_model(X):
    #train_size=int(len(X)*0.75)
    #train,test=X[0:train_size], X[train_size:]
    model=HWES(X,trend='add')
    fitted=model.fit(optimized=True,use_brute=True)
    pred=fitted.forecast(steps=100)
    return pred


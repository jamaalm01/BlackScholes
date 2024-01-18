#!/usr/bin/env python
# coding: utf-8

# In[21]:


#Implementation of the black-scholes formula

import numpy as np
from scipy.stats import norm

#defining variables

r = 0.01
S = 30
K = 40
T = 280/365
sigma = 0.30

#Calculates Black Scholes option price and greeks for a call or put
def blackScholes(r,S,K,T,sigma,callOrPut="c"):
    d1 = (np.log(S/K) + (r + sigma**2/2)*T)/(sigma*np.sqrt(T))
    d2 = d1 - sigma*np.sqrt(T)
    
    try:
        if callOrPut == "c":
            price = S*norm.cdf(d1,0,1) - K*np.exp(-r*T)*norm.cdf(d2,0,1)
            delta = norm.cdf(d1,0,1)
            gamma = norm.pdf(d1,0,1)/(S*sigma*np.sqrt(T))
            vega = 0.01*S*norm.pdf(d1,0,1)*np.sqrt(T)
            theta = (-S*norm.pdf(d1,0,1)*sigma/(2*np.sqrt(T)) - r*K*np.exp(-r*T)*norm.cdf(d2,0,1))/365
            rho = 0.01*K*T*np.exp(-r*T)*norm.cdf(d2,0,1)
        elif callOrPut == "p":
            price = K*np.exp(-r*T)*norm.cdf(-d2,0,1) - S*norm.cdf(-d1,0,1)
            delta = -norm.cdf(-d1,0,1)
            gamma = norm.pdf(d1,0,1)/(S*sigma*np.sqrt(T))
            vega = 0.01*S*norm.pdf(d1,0,1)*np.sqrt(T)
            theta = (-S*norm.pdf(d1,0,1)*sigma/(2*np.sqrt(T)) + r*K*np.exp(-r*T)*norm.cdf(-d2,0,1))/365
            rho = 0.01*-K*T*np.exp(-r*T)*norm.cdf(-d2,0,1)
        return print("Price:",round(price,2),"\nDelta:",round(delta,2),"\nGamma:",round(gamma,2),"\nVega:",round(vega,2),"\nTheta:",round(theta,3),"\nRho:",round(rho,2))
    except:
        print("Please confirm all option parameters above")


# In[22]:


blackScholes(r,S,K,T,sigma,"p")


#!/usr/bin/env python
# coding: utf-8

# In[3]:


import numpy as np
def FTBS1(initial_array,diffx,difft,windspeed):
    c=difft*windspeed/diffx
    q=np.zeros((5,5))
    q[0][0]=initial_array[0]
    q[0][1]=initial_array[1]
    q[0][2]=initial_array[2]
    q[0][3]=initial_array[3]
    q[0][4]=initial_array[4]
    print (q)
    q[1][1]=q[0][1]-c*(q[0][1]-q[0][0])
    q[1][2]=q[0][2]-c*(q[0][2]-q[0][1])
    q[1][3]=q[0][3]-c*(q[0][3]-q[0][2])
    q[1][4]=q[0][4]-c*(q[0][4]-q[0][3])
    q[2][2]=q[1][2]-c*(q[1][2]-q[1][1])
    q[2][3]=q[1][3]-c*(q[1][3]-q[1][2])
    q[2][4]=q[1][4]-c*(q[1][4]-q[1][3])
    
    return q
FTBS1([2,1,2,1,2],1,0.1,5)


# In[4]:


import numpy as np
def FTBS2(initial_array,diffx,difft,windspeed):
    c=difft*windspeed/diffx
    m=len(initial_array)
    q=np.zeros((m,m))
    for i in range (0,m):
        q[0][i]=initial_array[i]
    for k in range (1,m):
        for h in range (k,m):
            q[k][h]=q[k-1][h]-c*(q[k-1][h]-q[k-1][h-1])
    return q
FTBS2([2,4,7,3,6,8],1,0.1,5)


# In[8]:


#Import necessary modules.
import numpy as np
import matplotlib.pyplot as plt
#Create function to call in initial values for q and elements for Courant number.
def drw(initial_array,diffx,difft,windspeed):
#Calculate Courant number.
    c=difft*windspeed/diffx
#Work out how much data can be generated, based on amount in input.
    m=len(initial_array)
#Set up array of zeros to fill.
    q=np.zeros((m,m))
#Generate data based on FTBS scheme.
    for i in range (0,m):
        q[0][i]=initial_array[i]
    for k in range (1,m):
        for h in range (k,m):
            q[k][h]=q[k-1][h]-c*(q[k-1][h]-q[k-1][h-1])
#Place q values for each matching time and space position of matrix in an array.
    qx=q.diagonal()
#Generate x values for graph.
    x=np.arange(0,diffx*m,diffx)
#Plot model.
    plt.plot(x,qx)
    plt.xlabel('x')
    plt.ylabel('q(x)')
    plt.title('Advection Equation FTBS')
    plt.show()
    return
drw([2,4,7,3,6,8],1,0.1,5)
drw([20,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],1,0.1,5)


# In[24]:


def FTCS1(initial_array,diffx,difft,windspeed):
    c=difft*windspeed/diffx
    if c>1 or c<0:
        q="Your Courant number is outside of the bounds for stability."
    else:
        m=len(initial_array)
        q=np.zeros((m,m))
        for i in range (0,m):
            q[0][i]=initial_array[i]
        for k in range (1,m-1):
            for h in range (k,m-1):
                q[k][h]=q[k-1][h]-(c/2)*(q[k-1][h+1]-q[k-1][h-1])
    return q
FTCS1([3,4,5,6],7,9,30)


# In[25]:


FTCS1([3,4,5,6],1,0.1,3)


# In[12]:


def LW(initial_array,diffx,difft,windspeed):
    c=difft*windspeed/diffx
    m=len(initial_array)
    q=np.zeros((m,m))
    for i in range (0,m):
        q[0][i]=initial_array[i]
    for k in range (1,m-1):
        for h in range (k,m-1):
            q[k][h]=q[k-1][h]-(c/2)*(q[k-1][h+1]-q[k-1][h-1])+(c**2/2)*(q[k-1][h+1]-2*q[k-1][h]+q[k-1][h-1])
    return q
LW([2,0,0,0,0,0],1,0.1,5)


# In[ ]:





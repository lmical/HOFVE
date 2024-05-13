import numpy as np
import matplotlib.pyplot as plt
from mpmath import mp



#==============================================================
#Function to get the evaluation of Lagrange polynomials associated to certain nodes
#in certain points x
#In, particular the evaluations of l_k are given
#==============================================================
def lagrange_basis(nodes,x,k):
    """
    INPUT
     nodes, nodal points
     x,     vector where to evaluate the Lagrangian function l_k
     k,     index of the Lagrangian function to be evaluated
    OUTPUT
     y,     vector of evaluations of l_k in the abscissae in x
    """

    y=np.zeros(len(x))
    for ix, xi in enumerate(x):
        tmp=[(xi-nodes[j])/(nodes[k]-nodes[j])  for j in range(len(nodes)) if j!=k]
        y[ix]=np.prod(tmp)
    return y
#==============================================================
#
#
#
#==============================================================
#Function to get the derivative of Lagrange polynomials associated to certain nodes
#in certain points x
#In, particular the evaluations of d l_k are given
#==============================================================
# NB: Derivatives compute in multiple precision and then turned into float
#==============================================================
def lagrange_deriv(nodes,x,k):
    """
    INPUT
     nodes, nodal points
     x,     vector where to evaluate the derivative of the Lagrangian function l_k
     k,     index of the Lagrangian function, whose derivative has to be evaluated
    OUTPUT
     y,     vector of evaluations of l_k in the abscissae in x
    """

    # Added to deal with P0
    if len(nodes)==1:
	    return np.zeros(len(x))
    
    y=np.zeros(len(x))
    for indi in range(len(x)):
        f=mp.mpf(0) 
        for j in range(len(nodes)):
            p=mp.mpf(1)
            if k!=j:
                for l in range(len(nodes)):
                    if l!=k and l!=j: 
                        p=p*(x[indi]-nodes[l])/(nodes[k]-nodes[l])
                f = f + p/(nodes[k]-nodes[j])
        f=float(f)
        y[indi]=f
    return y
#==============================================================
#
#
#
#==============================================================
def convert_vector_mp_2_np(v):
    """
    From arbitrary precision to numpy array
    """
    vnp = np.zeros(len(v))
    for i in range(len(v)):
        vnp[i]=float(v[i])
    return vnp

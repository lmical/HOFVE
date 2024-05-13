import numpy as np
import reference_element

#INPUTS
p  = 5   #Number of cells in the single low order stencil, as well as number of stencils
xi = 1 #Quadrature point in [0,1]

#Actual position
xi=xi+(p-1.0) #I need to shift

#Interfaces in the stencils for interpolation
interfaces_HO=np.arange(0,2*p)
interfaces_LO=np.zeros((p,p+1)) #Row stencil, Column interfaces
for inds in range(p):
    interfaces_LO[inds,:]=np.arange(inds,inds+p+1)

# print(interfaces_HO)
# print(len(interfaces_HO))
# print(interfaces_LO)
# print(interfaces_LO.shape)

#Values derivatives of lagrangian polynomials in xi
derivatives_HO=np.zeros(2*p)
for indk in range(2*p): #2p interfaces
    derivatives_HO[indk]=reference_element.lagrange_deriv(interfaces_HO,[xi],indk)[0]

derivatives_LO=np.zeros((p,p+1))
for inds in range(p):
    for indk in range(p+1):
        derivatives_LO[inds,indk]=reference_element.lagrange_deriv(interfaces_LO[inds,:],[xi],indk)[0]

# print(derivatives_HO)
# print(np.sum(derivatives_HO))
# print(derivatives_LO)
# for inds in range(p):
#     print(np.sum(derivatives_LO[inds,:]))

#Coefficients of u_j in the HO polynomial
coefficients_HO=np.zeros(2*p-1)
for indk in range(2*p-1):
    coefficients_HO[indk]=np.sum(derivatives_HO[indk+1:])
    # print(derivatives_HO[indk+1:])


#Coefficients of u_j in the LO polynomials
coefficients_LO=np.zeros((p,p))
for inds in range(p): #Loop on the stencils
    for indk in range(p): #Loop on the cells
        coefficients_LO[inds,indk]=np.sum(derivatives_LO[inds,indk+1:])
        # print(derivatives_LO[inds,indk+1:])


# print(coefficients_HO)
# print(coefficients_LO)

#Linear system
#Naively speaking we have that alpha*coefficients_LO=coefficients_HO

#Matrix of the coefficients of the linear system
M=np.zeros((2*p-1,p)) #Rows are the u_j #Columns are the stencils
#Rhight-hand side of the linear system
r=np.zeros(2*p-1) #Rows are the u_j

end=-1
for inds in range(p): #Loop on the stencil
    for indk in range(p): #Loop on the cells of the stencil
        #Global index of the cell
        indr=inds+indk
        M[indr,inds]=coefficients_LO[inds,indk]

r=coefficients_HO.copy()

x=np.linalg.lstsq(M,r,rcond=None)
print(x[0])
# print(M@x[0]-r)
# HOFVE 
High Order Finite Volume solver for the Euler equations

# Flags
SW           -> Runs in SW mode: p=K*ro^\gamma (equivalence for K=g/2 and \gamma=2)<br />
WELLBALANCED -> Subtracts the space residual of the known steady state to have WB with respect to it<br />
PATANKAR     -> modified Patankar on the density only

# Space discretizations
1              = FV first order<br />
2=20           = standard MUSCL<br />
21,22,23,24,25 = different MUSCL schemes (from the picture MUSCL_scehemes.png)<br />
3              = W3<br />
4              = W5

# Time schemes
1 digit: 1 explicit euler, 2 SSPRK2, 3 SSPRK3, 4 SSPRK64,  5 RK65<br />
2 digits: 1n DeCn, 2n mPDeCn (NB: mPDeCn to be run with PATANKAR flag)


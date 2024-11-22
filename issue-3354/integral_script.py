# here a quick dirty python implementation of the integral in appendix A
import numpy as np


def lambda_by_another_name(tau, omega):
    """
    words
    """
    p = 1.0 - omega**2.0

    if (p < 0):
        integral = (1.0 / np.sqrt(np.abs(p)))*np.arcsin((1.0+omega*tau)/(tau+omega))
    #elif(p==0): # add and almost for # but how does htisfnkjve vkeqf vf # i dont this can be
    #    integral =
    else:
        integral = (1.0 / np.sqrt(np.abs(p)))*np.log((2.0*(1.0+tau*omega-np.sqrt(p*(1-tau**2.0))))/(tau+omega))

    return integral



# input to the function
Ro = 13.09 #14.016 # outer radius of the toroidal structure
Ri = 4.45 #3.911 # inner radius of the toroidal structure
Rm = 7.88 # peak R posiiton of the toroidal structre
H = 7.90 #8.089 #Max height of toroidal structure

#question is how do i select theta1?


theta1 = np.pi * (1.0/180.0)
theta2 = np.pi / 2.0 + theta1
a = (Ro - Ri) / 2.0
Rbar = (Ro + Ri) /2.0
A = Rbar / a
delta = (Rbar - Rm ) /a
kappa = H / a
iota = (1.0 + delta) / kappa

denom = np.cos(theta1) + np.sin(theta1) - 1.0

R1 = H*((np.cos(theta1) + iota*(np.sin(theta1) - 1.0)) / denom)
R2 = H*((np.cos(theta1) - 1.0 + iota* np.sin(theta1)) / denom)
R3 = H * (1 - delta) / kappa

Rc1 = (H / kappa) * (A  + 1.0) - R1
Rc2 = Rc1 + (R1 - R2) * np.cos(theta1)
Rc3 = Rc2
Zc1 = 0.0
Zc2 = (R1 - R2) * np.sin(theta1)
Zc3 = Zc2 + R2 - R3
# i assume index 4,5,6 are reflexed around z = 0  wee need to add these now
Rc4 = Rc1
Rc5 = Rc2
Rc6 = Rc3
Zc4 = -Zc1
Zc5 = -Zc2
Zc6 = -Zc3

#tau = np.array([[np.cos(theta1),np.cos(theta1+theta2),-1.0,np.cos(theta1),np.cos(theta1+theta2),-1.0],[1.0,np.cos(theta1),np.cos(theta1+theta2),1.0,np.cos(theta1),np.cos(theta1+theta2)]]) # this does habe k up to 6 indices
tau = np.array([[np.cos(theta1),np.cos(theta1+theta2),-1.0],[1.0,np.cos(theta1),np.cos(theta1+theta2)]]) # this does habe k up to 6 indices

#omega = np.array([Rc1/R1,Rc2/R2,Rc3/R3,Rc4/R1,Rc5/R2,Rc6/R3])
omega = np.array([Rc1/R1,Rc2/R2,Rc3/R3])

print(tau.shape)
chi1 = (Zc3 + np.abs(-Zc3)) / Ri
chi2 = 0.0
print(chi1)
for k in range(len(omega)):
    chi2 = chi2 + np.abs(lambda_by_another_name(tau[1,k],omega[k])-lambda_by_another_name(tau[0,k],omega[k]))
    print("chi2 =", chi2)
    print("lambda_1,",k," =", lambda_by_another_name(tau[1,k],omega[k]))
    print("lambda_0,",k," =", lambda_by_another_name(tau[0,k],omega[k]))

print(chi1+chi2)
phi = (chi1 + 2.0 * chi2) / (2.0 * np.pi)

print(phi)



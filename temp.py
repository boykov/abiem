j = 0
print "sigma(1) = ", self.P.sigma[j]
print "phi = ", self.P.intphi_over[j]
limA = self.P.intphi_over[j]/((2*math.pi)**(1.5)*self.P.sigma[j])
print "phi*limA = ", limA
phi_r = self.P.intphi_under[j]/(4*math.pi)
print "phi/(4pi*r) = ", phi_r
print "*/* =", limA/phi_r
print "sum = ", sum(self.P.intphi_under[:])/(4*math.pi)
print "φ(r)/r full =", 0.691797810029
print "φ(r)/r - sum =", sum(self.P.intphi_under[:])/(4*math.pi) - self.P.gauss[0,2]

self.centres[:] = self.quadphi_under[:,0]
self.C[:] = self.quadphi_under[:,1]
integ.calcomp2()

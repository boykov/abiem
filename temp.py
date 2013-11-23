j = 20
print "sigma(1) = ", self.P.sigma[j]
print "phi = ", self.P.intphi_over[j]
limA = self.P.intphi_over[j]/((2*math.pi)**(1.5)*self.P.sigma[j])
print "phi*limA = ", limA
phi_r = self.P.intphi_under[j]/(4*math.pi)
print "phi/(4pi*r) = ", phi_r
print "*/* =", limA/phi_r

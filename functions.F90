!- Mathematical functions used by endgame
module functions
 use config,    only:r0,T0,cs0,sigma_g0,p,q,mstar,rho_g0,eta0,rbump,epsi0,ibump,phi,w,epsimax,alpha
 use config,     only:gg,pi,au,earthr

 implicit none

 public         :: sigma_g,Temp,omega_k,cs,vdrift,vvisc,rho_g,vk,hoverr,h,epsi

 contains

real function sigma_g(r)
 real, intent(in) :: r

 sigma_g = sigma_g0*(r/r0)**(-p) + ibump*phi*sigma_g0*exp(-(r-rbump)**2/(2*w**2))
 return
end function sigma_g

real function rho_g(r)
 real, intent(in) :: r

 rho_g = sigma_g(r)/(sqrt(2*pi)*h(r))
 return
end function rho_g

real function Temp(r)
 real, intent(in) :: r

 Temp = T0*(r/r0)**(-q)
 return
end function Temp

real function omega_k(r)
 real, intent(in) :: r

 omega_k = sqrt(gg*mstar/r**3)
 return
end function omega_k

real function cs(r)
 real, intent(in) :: r

 cs = cs0*(r/r0)**(-q/2.)
 return
end function cs

real function vk(r)
 real, intent(in) :: r

 vk = omega_k(r)*r
 return
end function vk

real function vdrift(St,r)
 real, intent(in) :: St,r

 vdrift = r/press(r)*hoverr(r)**2*St/((1+epsi(r))**2+St**2)*vk(r)*&
(press(r+earthr)-press(r-earthr))/(2*earthr)
 return
end function vdrift

real function vvisc(St,r)
 real, intent(in) :: St,r
 real             :: denom,num,tmp1,tmp2 
 denom = 0.5*r*rho_g(r)*vk(r)
 tmp1 = -1.5*rho_g(r-earthr)*nu(r-earthr)*(r-earthr)**2*omega_k(r-earthr)
 tmp2 = -1.5*rho_g(r+earthr)*nu(r+earthr)*(r+earthr)**2*omega_k(r+earthr)
 num = (tmp2-tmp1)/(2*earthr)
 vvisc = num/denom * (1+epsi(r))/((1+epsi(r))**2 + St**2)
 return
end function vvisc

real function nu(r)
 real, intent(in) :: r
 nu = alpha*cs(r)*h(r)
 return
end function nu

real function hoverr(r)
 real, intent(in) :: r

 hoverr = h(r)/r
 return
end function hoverr

real function h(r)
 real, intent(in) :: r

 h = cs(r)/omega_k(r)
 return
end function h

real function press(r)
 real, intent(in) :: r

 press = cs(r)**2*rho_g(r)
 return
end function press

real function epsi(r)
 real, intent(in) :: r

 epsi = epsi0 + ibump*(epsimax-epsi0)*exp(-(r-rbump)**2/(2*w**2))
 return
end function epsi

end module functions

!- Mathematical functions used by endgame
module functions
 use config,    only:r0,T0,cs0,sigma_g0,p,q,mstar,rho_g0,eta0,rsnow,epsi0,ibump,a,b,epsimax
 use config,     only:gg,pi,au,earthr

 implicit none

 public         :: sigma_g,Temp,omega_k,cs,vdrift,rho_g,vk,hoverr,h,epsi

 contains

real function sigma_g(r)
 real, intent(in) :: r

 sigma_g = sigma_g0*(r/r0)**(-p) + ibump*sigma_g0/a*exp(-b*((r-rsnow)/rsnow)**2)
 return
end function sigma_g

real function rho_g(r)
 real, intent(in) :: r

 rho_g = sigma_g(r)/(2*pi*h(r))
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

 epsi = (epsi0 + ibump*(epsimax-epsi0)*exp(-b*((r-rsnow)/rsnow)**2))
 return
end function epsi

end module functions

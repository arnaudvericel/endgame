module evolve
 implicit none
 contains

subroutine evol(r,s,dsdt,vd,vdri,vvi,St,vrelonvfrag,rho,iam,iwas)
 use config,     only:dt,alpha,ifrag,vfrag,smin,vfragin,vfragout,&
                      rsnow,Tsnow,isnow,istate,igrow,Ro,racc
 use functions,  only:rho_g,cs,omega_k,vdrift,vvisc,Temp,epsi!,vset

 real, intent(inout)    :: r,s,rho
 real, intent(out)      :: St,dsdt,vrelonvfrag,vd,vdri,vvi

 integer, intent(inout) :: iam,iwas

 real                   :: vrel,ts

 ts = rho*s / (rho_g(r) * cs(r))
 St = ts * omega_k(r)

 if (igrow == 1) then
    vrel = sqrt(2**(3./2.)*alpha*Ro) * cs(r) * sqrt(St)/(1+St)
    dsdt = epsi(r)*rho_g(r)/rho*vrel
    if (ifrag == 0) then
       vrelonvfrag = vrel
    else
       select case(isnow)
       case(0)
           vrelonvfrag = vrel/vfrag
       case(1)
           if (r<=rsnow) then
              vrelonvfrag = vrel / vfragin
              iam = 1
           else
              vrelonvfrag = vrel / vfragout
              iam = 2
           endif
       case(2)
           if (Temp(r) >= Tsnow) then
              vrelonvfrag = vrel / vfragin
              iam = 1
           else
              vrelonvfrag = vrel / vfragout
              iam = 2
           endif
       end select
    endif
    if (vrelonvfrag < 1. .or. ifrag == 0) s = s + dsdt*dt
    if (vrelonvfrag >= 1. .and. ifrag > 0) s = s - dsdt*dt
    if (s < smin) s = smin
 else
    vrelonvfrag = 0.
    dsdt        = 0.
 endif

 if (istate == 1) then
    if (iam - iwas == -1) call sublimate(s,rho)
    if (iam - iwas == 1) call condense(s,rho,r)
    iwas = iam
 endif

 vdri = vdrift(St,r)
 vvi  = vvisc(St,r)
 vd = vdri + vvi
 if (r > racc) r  = r + vd*dt

end subroutine evol

subroutine condense(s,rho,r)
 use config,     only:rho1,rho2,abun
 use functions,  only:epsi
 real,intent(inout)    :: s,rho
 real,intent(in)       :: r

 real                  :: lambda

 lambda = abun/epsi(r)*rho1/rho2

 s = s * (1 + lambda)**(1./3.)
 rho = (rho + lambda*rho2) / (1 + lambda)

end subroutine condense

subroutine sublimate(s,rho)
 use config,     only:rho1,rho2,mratio,smin
 real,intent(inout)    :: s,rho

 s = s * (rho2*(1-mratio) / (rho1*mratio + rho2*(1-mratio)))**(1./3.)
 if (s < smin) s = smin
 rho = rho1

end subroutine sublimate

end module evolve

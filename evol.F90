   module evolve
implicit none
contains
!-- initialise run
subroutine init()
 use config,     only:p,q,mdisc,mstar,rin,rout,r0,T0,mu,epsi0,rho1,rho2,ifrag,isnow,rsnow,&
Tsnow,nsteps,nmax,sigma_g0,cs0,dt,tmax,rho_g0,vfrag,vfragin,&
vfragout,alpha,smin,iam,iwas,rho,abun,mratio,ndust,output,iamhere,ndumps,istate,&
rho,s,r,ibump,phi,w,epsimax,iam,iwas,igrow,disc,toto,dir
 use config,      only:pi,au,solarm,kboltz,mh
 use functions,  only:omega_k,vk,cs,h,hoverr,press,Temp,epsimax

 integer             :: i = 0
 logical             :: iexist = .false.

 !- Read from disc.in file
 inquire(file='disc.in', exist=iamhere)
 if (iamhere .eqv. .false.) print*, 'Error: input file disc.in not found'
 open(unit=1,file="disc.in",form="formatted",status="old",action="read")
 read(1,*)
 read(1,*)
 read(1,*) p,q,mdisc,mstar,rin,rout,r0,T0,mu,alpha,phi,w
 read(1,*)
 read(1,*)
 read(1,*)
 read(1,*) epsi0,rho1,rho2,smin,epsimax
 read(1,*)
 read(1,*)
 read(1,*)
 read(1,*) vfrag,vfragin,vfragout,rsnow,Tsnow,abun,mratio
 read(1,*)
 read(1,*)
 read(1,*)
 read(1,*) igrow,ifrag,isnow,ibump,istate
 read(1,*)
 read(1,*)
 read(1,*)
 read(1,*) nsteps,nmax,ndumps
 read(1,*)
 read(1,*)
 read(1,*)
 read(1,*) disc
 close(1)

 inquire(file=disc, exist=iexist)
 
 if (.not. iexist) then
    call system('mkdir '// adjustl(trim(disc)))
    write(*,*) "Creating directory ",adjustl(trim(disc))
 else
    write(*,*) "Directory ",adjustl(trim(disc)), " already present, skipping creation"
 endif

 write(toto(1),'(a,i1)') 'g',igrow
 if (ifrag == 1) then
    if (isnow == 0) write(toto(2),'(a,i1,a,i2,a)') 'f',ifrag,'V',int(vfrag),' '
    if (isnow == 1) then
       write(toto(2),'(a)') 'f1'
       write(toto(3),'(a,i1,a,i3,a,i2,a,i2)') 's',isnow,'r',int(rsnow),'V',int(vfragin),'-',int(vfragout)
    endif
    if (isnow == 2) then
       write(toto(3),'(a,i1,a,i3,a,i2,a,i2)') 's',isnow,'T',int(Tsnow),'V',int(vfragin),'-',int(vfragout)
       write(toto(2),'(a)') 'f1'
    endif
 else
    write(toto(2),'(a)') ' '
    write(toto(3),'(a)') ' '
 endif
 if (ibump == 0) write(toto(4),'(a)') ' '
 if (ibump == 1) write(toto(4),'(a,i1,a,f2.1,a,i4)') 'b',ibump,'phi',phi,'w',int(w)
 if (istate == 0) write(toto(5),'(a)') ' '
 if (istate == 1) write(toto(5),'(a,i1,a,f3.2)') 'st',istate,'mr',mratio

 do i=1,5
    call stripspaces(toto(i))
 enddo

 write(dir,*) adjustl(trim(disc)) // '/' // trim(toto(1)) // trim(toto(2)) // trim(toto(3)) // trim(toto(4)) // trim(toto(5))

 inquire(file=adjustl(trim(dir)), exist=iexist)
 if (.not. iexist) then
    call system('mkdir ' // dir)
    write(*,*) "Creating directory ",dir
 else
    write(*,*) "Directory ",adjustl(trim(dir)), " already present, skipping creation"
 endif
 !- Initialise parameters using the infile
 mdisc = mdisc * solarm
 mstar = mstar * solarm
 rin = rin * au
 rout = rout * au
 r0 = r0 * au
 w = w * au
 rsnow = rsnow * au
 sigma_g0 = mdisc/(2*pi) / (r0**2/(2-p)*((rout/r0)**(2-p) &
- (rin/r0)**(2-p)) + ibump*sqrt(pi)/4*(erf((rout-rsnow)/rsnow)-erf((rin-rsnow)/rsnow)))
 cs0 = sqrt(kboltz*T0/(mu*mh))
 rho_g0 = sigma_g0 / (sqrt(2*pi)*cs0/omega_k(r0))
 rho1 = rho1 * 1000
 rho2 = rho2 * 1000
 dt = 2*pi/omega_k(r0) / nsteps
 tmax = dt * nsteps * nmax
 iam(:) = 0
 iwas(:) = 0

 !- Read from dust.in file
 inquire(file='dust.in', exist=iamhere)
 if (iamhere .eqv. .false.) print*, 'Error: input file dust.in not found'
 open(unit=70,file="dust.in",form="formatted",status="old",action="read")
 read(70,*)
 read(70,*) ndust
 read(70,*)
 read(70,*)
 do i=1,ndust
    read(70,*) s(i),r(i),rho(i)
    r(i) = r(i) * au
    if (isnow > 0) rho(i) = (mratio*rho1 + (1-mratio)*rho2)
    if (isnow == 0) rho(i) = rho(i) * 1000
    iam(i) = 0
    iwas(i) = 0
 enddo
 close(80)

!- initialise dump files
 do i=1,ndust
    if (i<10) then
       write(output(i),'(a,i1,a)') 'p00',i,'.dat'
    elseif (i<100) then
       write(output(i),'(a,i2,a)') 'p0',i,'.dat'
    else
       write(output(i),'(a,i3,a)') 'p',i,'.dat'
    endif
 enddo

 !- write dust data file
 inquire(file='dust.dat', exist=iamhere)
 if (iamhere) open(unit=666,file="dust.dat",form="formatted",status="old",action="write")
 if (.not. iamhere) open(unit=666,file="dust.dat",form="formatted",status="new",action="write")
 do i=1,ndust
    write(666,'(a10,es10.2,i6,f10.2)') output(i),s(i),int(r(i)/au),rho(i)
 enddo
 close(666)

 write(6,'(a,i3,a)')'The Physical properties of your disc at ',int(r0/au), ' AU:'
 write(6,'(a,es8.2,a)')'sigma_gas = ',sigma_g0, ' kg/m2'
 write(6,'(a,es8.2,a)')'cs        = ',cs0,      ' m/s'
 write(6,'(a,es8.2,a)')'rho_gas   = ',rho_g0,   ' kg/m3'
 write(6,'(a,es8.2,a)')'P_gas     = ',press(r0),' Pa'
 write(6,'(a,es8.2,a)')'H/R       = ',hoverr(r0)

 write(6,'(a,i3,a)')'Currently ',ndust, ' grain(s) in the disc.'

#ifdef THANOS
 write(6,*) 'Your disc is perfectly balanced, as all things should be.'
#endif

end subroutine init

subroutine evol(r,s,dsdt,vd,vdri,vvi,St,vrelonvfrag,rho,iam,iwas)
 use config,     only:dt,alpha,ifrag,vfrag,smin,vfragin,vfragout,&
rsnow,Tsnow,isnow,istate,igrow
 use functions,  only:rho_g,cs,omega_k,vdrift,vvisc,Temp,epsi!,vset
 use config,     only:Ro

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
 r  = r + vd*dt

end subroutine evol

subroutine condense(s,rho,r)
 use config,     only:rho1,rho2,abun,smin
 use functions,  only:epsi
 real,intent(inout)    :: s,rho
 real,intent(in)       :: r

 real                  :: lambda

 lambda = abun/epsi(r)*rho1/rho2

 s = s * (1 + lambda)**(1./3.)
 if (s < smin) s = smin
 rho = (rho + lambda*rho2) / (1 + lambda)

end subroutine condense

subroutine sublimate(s,rho)
 use config,     only:rho1,rho2,mratio,smin
 real,intent(inout)    :: s,rho

 s = s * (mratio*rho2 / ((1-mratio)*rho1 + mratio*rho2))**(1./3.)
 if (s < smin) s = smin
 rho = rho1

end subroutine sublimate

subroutine stripspaces(randomstring)
 character(len=*),intent(inout) :: randomstring
 integer                        :: stringlen
 integer                        :: last,actual

 stringlen = len(randomstring)
 last      = 1
 actual    = 1

 do while (actual < stringlen)
    if (randomstring(last:last) == ' ') then
       actual = actual + 1
       randomstring(last:last) = randomstring(actual:actual)
       randomstring(actual:actual) = ' '
    else
       last = last + 1
       if (actual < last) actual = last
    endif
end do

end subroutine

end module evolve

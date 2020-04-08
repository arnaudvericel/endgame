module initial
 implicit none
 contains
 
 subroutine init()
 use config
 use functions,  only:omega_k,vk,cs,h,hoverr,press,Temp,epsimax

 integer             :: i = 0
 logical             :: exist = .false.
 real                :: titeuf

 !- Read from disc.in file
 inquire(file='disc.in', exist=iamhere)
 write(*,*) "Reading input file disc.in"
 if (iamhere .eqv. .false.) print*, 'Error: input file disc.in not found'
 open(unit=1,file="disc.in",form="formatted",status="old",action="read")
 read(1,*)
 read(1,*)
 read(1,*) p,q,mdisc,mstar,racc,rin,rout,r0,T0,mu,alpha,phi,w
 read(1,*)
 read(1,*)
 read(1,*)
 read(1,*) epsi0,rho1,rho2,smin,epsimax
 read(1,*)
 read(1,*)
 read(1,*)
 read(1,*) vfrag,vfragin,vfragout,rsnow,rbump,Tsnow,abun,mratio
 read(1,*)
 read(1,*)
 read(1,*)
 read(1,*) igrow,ifrag,isnow,ibump,istate,ibr
 read(1,*)
 read(1,*)
 read(1,*)
 read(1,*) nsteps,nmax,ndumps
 read(1,*)
 read(1,*)
 read(1,*)
 read(1,*) disc, isort
 close(1)

 
 if (isort == 1) then
    inquire(file=disc, exist=exist)
    if (.not. exist) then
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

    inquire(file=adjustl(trim(dir)), exist=exist)
    if (.not. exist) then
       call system('mkdir ' // dir)
       write(*,*) "Creating directory",dir
    else
       write(*,*) "Directory ",adjustl(trim(dir)), " already present, skipping creation"
    endif
 else
    inquire(file="buffer", exist=exist)
    if (.not. exist) then
       call system('mkdir buffer')
       write(*,*) "Creating directory buffer"
       dir = "buffer/"
    else
       write(*,*) "Buffer already present"
    endif
 endif
 
 !- Initialise parameters using the infile
 mdisc = mdisc * solarm
 mstar = mstar * solarm
 rin = rin * au
 rout = rout * au
 r0 = r0 * au
 w = w * au
 rsnow = rsnow * au
 rbump = rbump * au
 if (.not. (abs(p-2.) < 1.e-5)) then
    sigma_g0 = mdisc/(2*pi) / (r0**2/(2-p)*((rout/r0)**(2-p) &
    - (rin/r0)**(2-p)) + ibump*sqrt(pi/2)*w*(erf((rbump-rin)/(sqrt(2.)*w))-erf((rbump-rout)/(sqrt(2.)*w))))
 else
    sigma_g0 = mdisc/(2*pi*r0**2*log(rout/rin) + ibump*sqrt(pi/2)*w*(erf((rbump-rin)/(sqrt(2.)*w))-erf((rbump-rout)/(sqrt(2.)*w))))
 endif
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
 write(*,*) "Reading input file dust.in"
 if (iamhere .eqv. .false.) print*, 'Error: input file dust.in not found'
 open(unit=70,file="dust.in",form="formatted",status="old",action="read")
 read(70,*)
 read(70,*) ndust
 read(70,*)
 read(70,*)
 do i=1,ndust
    read(70,*) titeuf,s(i),r(i),rho(i)
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

 write(6,'(a,i3,a)')'The Physical properties of the disc at ',int(r0/au), ' AU:'
 write(6,'(a,es8.2,a)')'sigma_gas = ',sigma_g0, ' kg/m2'
 write(6,'(a,es8.2,a)')'cs        = ',cs0,      ' m/s'
 write(6,'(a,es8.2,a)')'rho_gas   = ',rho_g0,   ' kg/m3'
 write(6,'(a,es8.2,a)')'P_gas     = ',press(r0),' Pa'
 write(6,'(a,es8.2,a)')'H/R       = ',hoverr(r0)

 write(6,'(a,i3,a)')'Currently ',ndust, ' grain(s) in the disc.'

#ifdef THANOS
 write(6,*) 'The disc is perfectly balanced, as all things should be.'
#endif

end subroutine init

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

end module initial
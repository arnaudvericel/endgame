!- module defining all needed variables and parameters
module config
 implicit none

 !- Disc related
 real,public              :: p
 real,public              :: q
 real,public              :: r0
 real,public              :: rin
 real,public              :: rout
 real,public              :: mdisc
 real,public              :: mstar
 real,public              :: mu
 real,public              :: alpha
 real,public              :: sigma_g0
 real,public              :: rho_g0
 real,public              :: T0
 real,public              :: cs0
 real,public              :: eta0
 real,public              :: a
 real,public              :: b

!- Dust related
 real,public              :: epsi0
 real,public              :: rho1
 real,public              :: rho2
 real,public              :: Tsnow
 real,public              :: rsnow
 real,public              :: vfrag
 real,public              :: vfragin
 real,public              :: vfragout
 real,public              :: smin
 real,public              :: abun
 real,public              :: mratio
 real,public              :: epsimax

 integer,public           :: ndust
 integer,public           :: iam(200)
 integer,public           :: iwas(200)

!- Timestep related
 real,public              :: dt
 real,public              :: tmax

 integer,public           :: nsteps
 integer,public           :: nmax
 integer,public           :: ndumps
 integer,public           :: step = 0
 integer,public           :: nprev = 0
 integer,public           :: ntic = 0


!- options control related
 integer,public           :: igrow
 integer,public           :: ifrag
 integer,public           :: isnow
 integer,public           :: istate
 integer,public           :: ibump

!- files related
 character(len=10),public :: output(200)
 character(len=30),public :: disc
 character(len=15),public :: toto(5)
 character(len=100),public :: dir

 logical,public           :: iexist(200)
 logical,public           :: iamhere

 integer,public           :: idir

!- actual variables that evolve during the simulation
 real,public              :: vd(200)
 real,public              :: vz(200)
 real,public              :: St(200)
 real,public              :: vrelonvfrag(200)
 real,public              :: r(200)
 real,public              :: z(200)
 real,public              :: s(200)
 real,public              :: dsdt(200)
 real,public              :: rho(200)

!**************
!  parameters
!**************

!- Mathematical constants
 real(kind=8), parameter  :: pi         =  3.1415926536d0

!- Physical constants
 real(kind=8), parameter  :: gg         = 6.672041d-11
 real(kind=8), parameter  :: kboltz     = 1.38066d-23
 real(kind=8), parameter  :: avogadro   = 6.0221408577d23
 real(kind=8), parameter  :: Ro         = 3.d0
 real(kind=8), parameter  :: mh         = 1.67e-27

!- Mass, radius and distance conversions
 real(kind=8), parameter  :: solarm     = 1.9891d30
 real(kind=8), parameter  :: solarr     = 6.959500d8
 real(kind=8), parameter  :: earthm     = 5.979d24
 real(kind=8), parameter  :: earthr     = 6.371315d6
 real(kind=8), parameter  :: au         = 1.496d11
 real(kind=8), parameter  :: km         = 1.d3

!- Time
 real(kind=8), parameter  :: seconds = 1.d0
 real(kind=8), parameter  :: minutes = 6.0d1
 real(kind=8), parameter  :: hours = 3.6d3
 real(kind=8), parameter  :: days = 8.64d4
 real(kind=8), parameter  :: years = 3.1556926d7
end module config

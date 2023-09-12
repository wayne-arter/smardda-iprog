module const_numphys_h

  use const_kind_m

  implicit none
  private


  real(kr4) , parameter, public :: const_pi= 3.14159265_kr4 !<   4 byte pi

  real(kr8) , parameter, public :: const_pid = 3.14159265358979_kr8 !<   8 byte pi

!> neutron mass etc
!! MCNP values
  real(kr8) , parameter, public :: &  !< local variable
 &const_aneut  = 1.008664967_kr8,            & !< Neutron mass in a.m.u..
 &const_avogad = 6.022043446928244e+23_kr8,     & !< Avogadro's number.
 &const_mneut  = .001_kr8*const_aneut/const_avogad,       & !< Neutron mass in kg
 &const_avgdn  = 1.0e-24_kr8*const_avogad/const_aneut,      & !< 1.e-24*Avogadro's number/neutron mass.
 &const_euler  = .577215664901532861_kr8,    & !< Euler constant used in elect ron transport.
 &const_fscon  = 137.0393_kr8,               & !< Inverse fine-structure const ant.
 &const_planck = 4.135732e-13_kr8,              & !< Planck constant.
 &const_rmass = 939.56563_kr8,              & !< rest mass of neutron in MeV
 &const_charge = 1.602176487e-19_kr8,     & !< absolute electron charge in C.
 &const_massu = 1.673774424e-27_kr8,        & !<  Hydrogen atom mass in kg.
 &const_massp = 1.6726219e-27_kr8,     & !< proton mass in kg.
 &const_masse = 9.10938e-31_kr8,           & !< mass of electron in kg.
 &const_slite  = 299.7925_kr8,                 & !< Speed of light in Mm/s (NB factor of 1e6)
 &const_epsilon0  = 8.8542e-12_kr8,                 & !< Permeability of free space
 &const_mu0  = 4*const_pid*1.e-7_kr8,                 & !< Permittivity of free space
 &const_kb = 1.38064852e-23_kr8      !< boltzmann constant

  real(kr4) , parameter, public :: const_pushinf = 1.e+8_kr4 !<  infinity for particle push

  real(kr4) , parameter, public :: const_pusheps = 1.e-8_kr4 !<  epsilon for particle push
  real(kr8) , parameter, public :: const_epsbdry = 1.e-6_kr8 !<  epsilon for particle push

  real(kr8) , parameter, public :: const_golden = (1+sqrt(5._kr8))/2 !<  golden ratio

  real(kr8) , parameter, public :: const_infty = 1.e+30_kr8 !<  numerical infinity
  real(kr4) , parameter, public :: const_degrad = const_pi/180 !<  degrees to radians conversion factor
  real(kr8) , parameter, public :: const_r36th = 1.0_kr8/36.0_kr8 !<  spl2d normaliser

end module const_numphys_h

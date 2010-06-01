module turbstats_vars
  implicit none
  ! 
  ! Meta info about the stats
  ! 

  integer, parameter :: nstats = 21
  ! column headings for table
  character(len=4 ), dimension(nstats) :: statid = (/&
       & 'Mn  ',&
       & 'Mi  ',&
       & 'Xcom',&
       & 'Ycom',&
       & 'Zcom',&
       & 'Vnx ',&
       & 'Vny ',&
       & 'Vnz ',&
       & 'Bx  ',&
       & 'By  ',&
       & 'Bz  ',&
       & 'Bix ',&
       & 'Biy ',&
       & 'Biz ',&
       & 'Bnx ',&
       & 'Bny ',&
       & 'Bnz ',&
       & 'Brms',&
       & 'Beta',&
       & 'Be_i',&
       & 'Be_n' &
       &/)
  !
  ! Ths stats themselves
  !
  real :: neutral_mass, ionized_mass
  ! Center of mass
  real :: xcom, ycom, zcom
  ! Globule neutral COM velocities
  real :: vnx, vny, vnz
  ! Mean field
  real :: bmeanx, bmeany, bmeanz
  real :: bmeanx_i, bmeany_i, bmeanz_i
  real :: bmeanx_n, bmeany_n, bmeanz_n
  real :: brms
  ! Position of ionization front
  real :: xifront
  ! velocity of neutral shell
  real :: vxshell
  ! average betas
  real :: beta, beta_i, beta_n


  namelist / allstats / &
       & neutral_mass, ionized_mass, &
       & xcom, ycom, zcom, vnx, vny, vnz,&
       & bmeanx, bmeany, bmeanz,&
       & bmeanx_i, bmeany_i, bmeanz_i,&
       & bmeanx_n, bmeany_n, bmeanz_n,&
       & brms,&
       & beta, beta_i, beta_n
end module turbstats_vars

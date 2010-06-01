module allstats_namelist
  implicit none
  namelist / allstats / &
       & neutral_mass, ionized_mass, &
       & xcom, ycom, zcom, vnx, vny, vnz,&
       & bmeanx, bmeany, bmeanz,&
       & bmeanx_i, bmeany_i, bmeanz_i,&
       & bmeanx_n, bmeany_n, bmeanz_n,&
       & brms, xifront
end module allstats_namelist

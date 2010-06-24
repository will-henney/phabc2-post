program globcollate
  use allstats_vars, only: nstats
  implicit none
  character(len=128) :: prefix
  integer, parameter :: maxtime = 400
  integer :: itime, krow = 0, j
  real, dimension(nstats+1,maxtime) :: table

  print *, 'File prefix?'
  read '(a)', prefix

  do itime = 1, maxtime
     call read_stats(itime)
  end do
  
  call write_table


contains

  subroutine write_table
    use allstats_vars, only: statid
    character(len=1) :: TAB = achar(9)
    character(len=256) :: hdrformat, dataformat
    integer :: k
    write(hdrformat, '(a,i0,3a)') '(', nstats+1, '(a,"', TAB, '"))'  
    write(dataformat, '(a,i0,3a)') '(', nstats+1, '(es10.3,"', TAB, '"))'  
    print *, trim(hdrformat)
    print *, trim(dataformat)
    open(12, &
         & file=trim(prefix)//'-'//'globstats-table.dat', &
         & action='write')
    write(12, trim(hdrformat)) '# Time', statid
    do k = 1, krow
       write(12, trim(dataformat)) table(:,k)
    end do
    close(12)
  end subroutine write_table
  
  subroutine read_stats(itime)
    use allstats_vars
    integer, intent(in) :: itime
    character(len=4) :: ctime
    integer :: istat
    write(ctime, '(i4.4)') itime
    istat = 0
    open(11, &
         & file=trim(prefix)//'-'//'globstats'//ctime//'.dat', &
         & action='read', status='old', iostat=istat)
    if (istat == 0) then
       read(11, allstats, iostat=istat)
       close(11)
    else
       ! Don't worry if we can't read the file - just skip this time
       return
    end if

    ! If we are here, that we have the stats for this time

    ! Increment the row counter
    krow = krow + 1


    ! we do it like this, so we don't need to worry about the column
    ! numbers, just about the order
    j = 1
    call addvar2row(real(itime))     ! first column is the time
    call addvar2row(neutral_mass)
    call addvar2row(ionized_mass)
    call addvar2row(xcom)
    call addvar2row(ycom)
    call addvar2row(zcom)
    call addvar2row(vnx)
    call addvar2row(vny)
    call addvar2row(vnz)
    call addvar2row(bmeanx)
    call addvar2row(bmeany)
    call addvar2row(bmeanz)
    call addvar2row(bmeanx_i)
    call addvar2row(bmeany_i)
    call addvar2row(bmeanz_i)
    call addvar2row(bmeanx_n)
    call addvar2row(bmeany_n)
    call addvar2row(bmeanz_n)
    call addvar2row(brms)
    call addvar2row(xifront)
    call addvar2row(vxshell)
    call addvar2row(beta)
    call addvar2row(beta_i)
    call addvar2row(beta_n)

  end subroutine read_stats

  subroutine addvar2row(var)
    real, intent(in) :: var
    table(j,krow) = var
    j = j + 1
  end subroutine addvar2row
      

end program globcollate

double precision function rtc()
  implicit none
  integer:: icnt,irate
  double precision, save:: scaling
  logical, save:: scale = .true.
  ! call CPU_TIME(rtc)
  call system_clock(icnt,irate)

  if(scale)then
     scaling=1.0/dble(irate)
     scale=.false.
  end if

  rtc = icnt * scaling

end function rtc

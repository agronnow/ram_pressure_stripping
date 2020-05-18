module sn_feedback_commons

  use amr_parameters
  use random

  implicit none

  integer,dimension(IRandNumSize) :: localseedsn=-1
  integer::iseed=0
#ifdef DELAYED_SN
  integer::nSN_prev = 0
  real(dp),dimension(NMAX_SN,3)::sn_coords
  real(dp),dimension(NMAX_SN)::rsn_sq
  integer,dimension(NMAX_SN)::sn_level
  integer,dimension(NMAX_SN)::sn_isrefined = 0

#endif

end module sn_feedback_commons

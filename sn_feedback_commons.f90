module sn_feedback_commons

  use amr_parameters
  use random

  implicit none

  integer,dimension(IRandNumSize) :: localseed=-1
  integer::iseed=0
#ifdef DELAYED_SN
  integer::nSN_prev = 0
  real(dp),dimension(NMAX_SN,3)::sn_coords
#endif

end module sn_feedback_commons

program diff_bim
use constants 
implicit none
integer, parameter :: dp = selected_real_kind(15,307)
! program to calculate ratediff and ratediss
real (dp) :: deltag,eta,ratediff,ratediss,temp,k
! read eta (in Pa*s=pascal*second) deltag (of reaction) and temperature (in K)
!read(*,*) eta,deltag,temp
read(*,*) eta,temp
ratediff=8*boltz*temp/3/eta*avog*1d3
k=exp(-deltag/r/temp)/ratm/temp
ratediss=ratediff*k
! rate
print*, ratediff
!print*, ratediss

end program diff_bim


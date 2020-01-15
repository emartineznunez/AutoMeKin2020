program termo
use constants
use atsymb
implicit none
integer, parameter:: b8 = selected_real_kind(14)
real(b8), parameter :: conv = 204548.28_b8
real, dimension (:),allocatable :: px,py,pz,w,vx,vy,vz
real :: desket,temp0,temp1,temp2
integer :: n,i,j,nate,ncycles,ii,step,delta
integer, dimension(:),allocatable :: inate
character*2,dimension(:),allocatable :: fasymb

! nate is the number of excited atoms
! inate are the indeces  of those atoms

read(*,*) delta
read(*,*) n,ncycles,step
read(*,*) nate
if(nate==0) then
   nate=n
   allocate(inate(nate))
   inate=(/ (i,i=1,n) /) 
else
   allocate(inate(nate))
   read(*,*) (inate(i),i=1,nate)
endif
read(*,*) n
allocate(px(n),py(n),pz(n),vx(n),vy(n),vz(n),w(n),fasymb(n))
!delta=ncycles*(step-1)
do ii=1+delta,ncycles+delta
   temp0=0.d0
   temp1=0.d0
   temp2=0.d0
   do i=1,n
      read(*,*) fasymb(i),vx(i),vy(i),vz(i)
      do j = 1 , natom
!         if(string_tolower(fasymb(i))==asymb(j)) w(i)=ams(j)
         if(fasymb(i)==asymb(j)) w(i)=ams(j)
      enddo
      px(i)=vx(i)*w(i)/conv
      py(i)=vy(i)*w(i)/conv
      pz(i)=vz(i)*w(i)/conv
      temp0=temp0+(px(i)**2+py(i)**2+pz(i)**2)/w(i)
      if(any(i .eq. inate)) then
        temp1=temp1+(px(i)**2+py(i)**2+pz(i)**2)/w(i)
      else
        temp2=temp2+(px(i)**2+py(i)**2+pz(i)**2)/w(i)
      endif
   enddo
   temp0=temp0/3/n/r
   temp1=temp1/3/nate/r
   if(nate<n) temp2=temp2/3/(n-nate)/r
   print*, ii,temp0,temp1,temp2
enddo

contains

function string_tolower( string ) result (new)
    character(len=*)           :: string

    character(len=len(string)) :: new

    integer                    :: i
    integer                    :: k

    new    = string
    do i = 1,len(string)
        k = iachar(string(i:i))
        if ( k >= iachar('A') .and. k <= iachar('Z') ) then
            k = k + iachar('a') - iachar('A')
            new(i:i) = achar(k)
        endif
    enddo
end function string_tolower


end program



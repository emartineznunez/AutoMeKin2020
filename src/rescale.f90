program termo
use atsymb
implicit none
integer, parameter:: b8 = selected_real_kind(14)
real(b8), parameter :: boltz = 0.00198717_b8
real(b8), parameter :: conv = 204548.28_b8
real, dimension (:),allocatable :: px,py,pz,w,vx,vy,vz
real :: desket,temp,temp0,factor,sumx,sumy,sumz,tempf,ener
integer :: n,i,j,nate
integer, dimension(:),allocatable :: inate,q
character*2,dimension(:),allocatable :: fasymb

! nate is the number of excited atoms
! inate are the indeces  of those atoms


read(*,*) n
allocate(q(3*n),px(n),py(n),pz(n),vx(n),vy(n),vz(n),w(n),fasymb(n))
do i=1,n
   read(*,*) fasymb(i),vx(i),vy(i),vz(i)
   do j = 1 , natom
!      if(string_tolower(fasymb(i))==asymb(j)) w(i)=ams(j)
      if(fasymb(i)==asymb(j)) w(i)=ams(j)
   enddo
enddo
read(*,*) temp
read(*,*) nate
if(nate==0) then
   nate=n
   allocate(inate(nate))
   inate=(/ (i,i=1,n) /)
else
   allocate(inate(nate))
   read(*,*) (inate(i),i=1,nate)
endif
ener=3*nate*0.5*boltz*temp
write(67,*) "TEmp=",temp
write(67,*) "Ener=",ener
write(67,*) "Number of atoms to be excited=",nate
write(67,*) "Atoms to be excited:",(inate(i),i=1,nate)

sumx = 0.d0
sumy = 0.d0
sumz = 0.d0
do i=1,n
   if(any(i .eq. inate) ) then
      px(i)=vx(i)*w(i)/conv
      py(i)=vy(i)*w(i)/conv
      pz(i)=vz(i)*w(i)/conv
      sumx = sumx + px(i)
      sumy = sumy + py(i)
      sumz = sumz + pz(i)
   else
      px(i)=0.
      py(i)=0.
      pz(i)=0.
   endif
enddo
sumx = sumx/nate
sumy = sumy/nate
sumz = sumz/nate

!    make sure that the total linear momentum equals zero
temp0=0.d0
do i=1,n
   if(any(i .eq. inate)) then
      px(i) = px(i) - sumx
      py(i) = py(i) - sumy
      pz(i) = pz(i) - sumz
   endif
   temp0=temp0+(px(i)**2+py(i)**2+pz(i)**2)/w(i)
end do
! temp0 is the initial temperature
temp0=temp0/3/nate/boltz
write(67,*) "Initial temp=",temp0
factor=sqrt(temp/temp0)

tempf=0.d0
ener=0.d0
do i=1,n
   if(any(i .eq. inate) ) then
      px(i) = px(i)*factor 
      py(i) = py(i)*factor 
      pz(i) = pz(i)*factor 
      tempf=tempf+(px(i)**2+py(i)**2+pz(i)**2)/w(i)
      ener=ener+0.5*(px(i)**2+py(i)**2+pz(i)**2)/w(i)
   endif
   vx(i)=px(i)/w(i)*conv
   vy(i)=py(i)/w(i)*conv
   vz(i)=pz(i)/w(i)*conv
   print*, vx(i),vy(i),vz(i)
end do
tempf=tempf/3/nate/boltz
write(67,*) "Final temp=",tempf
write(67,*) "Final ener=",ener,3*nate*0.5*boltz*tempf

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



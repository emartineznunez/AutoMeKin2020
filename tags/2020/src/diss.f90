program diss
use atsymb
implicit none
! program to rotate a molecule about its euler angles
! 1 is big
! 2 is small
integer :: i,n1,n2,n,j,iii,at1,at2
integer, dimension(:),allocatable :: iat1 
real :: xcm1,ycm1,zcm1,xcm2,ycm2,zcm2,wt,x,y,z
real :: dist,xdum,ydum,zdum,vx,vy,vz,delta,xref1,yref1,zref1,xref2,yref2,zref2

real,dimension(:),allocatable :: q,qq,w
character*2,dimension(:),allocatable :: fasymb
read(*,*) iii
read(*,*) n1,n2
allocate(iat1(n1))
read(*,*) (iat1(i),i=1,n1)
n=n1+n2
read(*,*) at1,at2
allocate(q(3*n),qq(3*n),fasymb(n),w(n))
do i=1,n
   read(*,*) fasymb(i),q(3*i-2),q(3*i-1),q(3*i)
   do j = 1 , natom
      if(fasymb(i)==asymb(j)) w(i)=ams(j)
   enddo
enddo
! com of second monomer (the small one)
xcm2=0
ycm2=0
zcm2=0
wt=0
do i=1,n
   if(all(i .ne. iat1)) then
      xcm2=xcm2+w(i)*q(3*i-2) 
      ycm2=ycm2+w(i)*q(3*i-1) 
      zcm2=zcm2+w(i)*q(3*i  ) 
      wt=wt+w(i)
   endif
enddo
xcm2=xcm2/wt
ycm2=ycm2/wt
zcm2=zcm2/wt
if(at2==-1) then
   xref2=xcm2
   yref2=ycm2
   zref2=zcm2
else 
   xref2=q(3*at2-2)
   yref2=q(3*at2-1)
   zref2=q(3*at2  )
endif
do i=1,n
   if(all(i .ne. iat1)) then
      qq(3*i-2)=q(3*i-2)-xref2
      qq(3*i-1)=q(3*i-1)-yref2
      qq(3*i  )=q(3*i  )-zref2
   endif
enddo

! com of first monomer (big one)
xcm1=0
ycm1=0
zcm1=0
wt=0
do i=1,n
   if(any(i .eq. iat1)) then
      xcm1=xcm1+w(i)*q(3*i-2) 
      ycm1=ycm1+w(i)*q(3*i-1) 
      zcm1=zcm1+w(i)*q(3*i  ) 
      wt=wt+w(i)
   endif
enddo
xcm1=xcm1/wt
ycm1=ycm1/wt
zcm1=zcm1/wt
if(at1==-1) then
   xref1=xcm1
   yref1=ycm1
   zref1=zcm1
else
   xref1=q(3*at1-2)
   yref1=q(3*at1-1)
   zref1=q(3*at1  )
endif
do i=1,n
   if(any(i .eq. iat1)) then
      qq(3*i-2)=q(3*i-2)-xref1
      qq(3*i-1)=q(3*i-1)-yref1
      qq(3*i  )=q(3*i  )-zref1
   endif
enddo

xdum=xref1-xref2
ydum=yref1-yref2
zdum=zref1-zref2
dist=sqrt(xdum*xdum + ydum*ydum + zdum*zdum)
vx=xdum/dist
vy=ydum/dist
vz=zdum/dist

delta=dist+(iii-1)*.2
print*, n
print*, delta,dist
DO I=1,n
   x=QQ(3*i-2)
   y=QQ(3*i-1)
   z=QQ(3*i  )
   if(any (i .eq. iat1)) then
      q(3*i-2)=x+delta*vx
      q(3*i-1)=y+delta*vy
      q(3*i  )=z+delta*vz
      print*, fasymb(i),q(3*i-2),q(3*i-1),q(3*i)
   else
      q(3*i-2)=x
      q(3*i-1)=y
      q(3*i  )=z
      print*, fasymb(i),q(3*i-2),q(3*i-1),q(3*i)
   endif
enddo

end program diss 

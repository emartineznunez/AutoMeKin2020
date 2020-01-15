program kmc
implicit none
! KMC computes the population vs time of a number of species involved in several dynamical processes
! rate: rate constant of a given processes  
! p:    population of a given species (p0 is its initial value)
! re (pr): reactant (product) for a given process
character*80 :: title
integer, parameter :: dp = selected_real_kind(15,307)
integer :: m,nran,nr,nesp,inran,i,j,k,l,mu,pd,kk,ia,iz,ijk,nespct
real(dp) :: t,tmax,tprint,tint,a0,r2a0,suma,vol,avog,conv
real(dp) :: rnd(2)
integer,dimension(:),allocatable :: re1,re2,pr1,pr2,n,iespct
real (dp) ,dimension(:),allocatable :: a,rate,c0,p,p0,cont
avog=6.022d23
read(*,"(a80)") title
print "(t3,a80)",title
read(*,*) m,nesp,nran
print*, m,nesp,nran
allocate(re1(m),re2(m),pr1(m),pr2(m),a(m),rate(m),p0(nesp),p(nesp),n(nesp),c0(nesp),cont(m))
n=(/ (l,l=1,nesp) /)
! The vol relates to population of species i p_i and concentration c_i as: vol=p_i/c_i/avog
! So, choose vol for a expected population of a given species i (p_i)
read(*,*) vol
conv=vol*avog
! nespct is the number of species with constant values of population throughout
! the simulation. For instance, in a catalyic cycles those species could be
! gases 
read(*,*) nespct
print*, "Number of species that keep their concentrations constant=",nespct
allocate(iespct(nespct))
do i=1,nespct
   read(*,*) iespct(i)
   print*, iespct(i)
enddo
print*, ""
print*, "Rates in (time-1),i.e.,s-1"
do i=1,m
   read(*,*) rate(i),re1(i),re2(i),pr1(i),pr2(i)
   if(re1(i)/=0.and.re2(i)/=0) rate(i)=rate(i)/conv
   print*, i,rate(i),re1(i),re2(i),pr1(i),pr2(i)
   if(re1(i)==re2(i)) rate(i)=2*rate(i)
enddo
print*, "Initial concentration (M) and population of each species"
do ijk=1,m
   cont(ijk)=0
enddo
do i=1,nesp
   read(*,*) c0(i)
   p0(i)=c0(i)*conv
   print*, i,c0(i),p0(i)
enddo
read(*,*) tmax,tint
ia=1
iz=nesp
print "(/,t3,a,1p,e10.2,a,/,t3,a,1p,e10.2,a,/)","Total time=",tmax," in input units","Step size =",tint," in input units"
big: do inran=1,nran
   print "(/,t3,a,i4,/,t3,a,500(i10))","Calculation number",inran, &
 "Time    ",(n(i),i=ia,iz)
   p=p0
   t=0.d0
   tprint=0.d0
   do while(tprint<tmax)
! keep the initial population for the iespct species
     do j=1,nespct
        p(iespct(j))=p0(iespct(j)) 
     enddo
     do j=1,m
        if(re1(j)==0) then
! unimolecular reaction
          a(j)=rate(j)*p(re2(j))
        else if(re2(j)==0) then
          a(j)=rate(j)*p(re1(j))
        else
          if(re1(j)/=re2(j)) a(j)=rate(j)*p(re1(j))*p(re2(j))
          if(re1(j)==re2(j)) a(j)=rate(j)*p(re1(j))*(p(re2(j))-1)/2
        endif
     enddo
     a0=sum(a)
     call random_number(rnd)
     t=t-log(rnd(1))/a0      
     do while (t>=tprint) 
        print "(e10.4,500(f15.0))",tprint,(p(i),i=ia,iz)
        tprint=tprint+tint
        if(tprint>tmax) cycle big
     enddo
     r2a0=rnd(2)*a0
     suma=0.d0
     s1: do mu=1,m
       suma=suma+a(mu)
       if(suma>=r2a0) exit s1
     enddo s1
     if(re1(mu)==0) then
! unimolecular reaction
        p(re2(mu))=p(re2(mu))-1
     else if(re2(mu)==0) then
        p(re1(mu))=p(re1(mu))-1
! unimolecular reaction
     else
        p(re1(mu))=p(re1(mu))-1
        p(re2(mu))=p(re2(mu))-1
     endif
     if(pr1(mu)==0) then
! only one product
        p(pr2(mu))=p(pr2(mu))+1
     else if(pr2(mu)==0) then
        p(pr1(mu))=p(pr1(mu))+1
     else 
        p(pr1(mu))=p(pr1(mu))+1
        p(pr2(mu))=p(pr2(mu))+1
     endif
     cont(mu)=cont(mu)+1  
   enddo
enddo big
print*,"Population of every species"
do i=1,nesp
   print "(i6,f20.0)",i,p(i)
enddo
print*,"counts per process"
do i=1,m
   print "(i6,f30.0)",i,cont(i)
enddo
end program kmc

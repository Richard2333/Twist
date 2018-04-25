!parameters from 10.1103/PhysRevB.88.125426
!needs layers file from overlapgen
!
module functions
   real(8),parameter :: d0=0.335
   real(8),parameter :: delt=0.0453
   real(8),parameter :: dela=0.142
   real(8),parameter :: t0=2.7,t12=-0.48
!~    real(8),parameter :: va1(3)=[0,1,0],va2(3)=[0,2,-1]
contains
   function dist(a,b) result(d)
      implicit none
      real(8) :: a(2),b(2),d
      d=abs(a(1)-b(1))+abs(a(2)-b(2))
   end function
   function hoppingvalue(v0,v1) result(t)
      implicit none
      real(8) :: t
      real(8) :: v0(3),v1(3),dv(3)
      real(8) :: d2,dz,dst
      dv=v1-v0
      d2=0.246*sqrt(dv(2)**2+dv(3)**2+dv(2)*dv(3))

      dz=0.335*dv(1)
      dst=sqrt(dz**2+d2**2)
      t=0
      if(abs(dv(1))<0.001) then
         t=t0*exp(-(dst-dela)/delt)
      else
         t= t12* exp(-(dst-d0)/delt)
      end if
      write(*,*) dv,dst
      return
   end function
end module
program main
   use functions
   implicit none
   integer,parameter :: nats=40

   real(8) :: atoms(nats,3),tm,tmc,j1(3),j2(3),va1(3),va2(3)
   integer :: i,j
   integer :: m,n
   integer :: cot
!~ 	real(8),allocatable:: hoppings
   i=1
   open(unit=234,file='layers.csv')
   open(unit=233,file='s1hr.dat')
   va1=0
   va2=0
   read(234,*) va1(2),va1(3)
   read(234,*) va2(2),va2(3)
   cot=0
   do i=1,nats
      read(234,*,end=1997) atoms(i,1),atoms(i,2),atoms(i,3)
      write (*,*) i
      cot=cot+1
   end do
!~ 100 continue
1997   write(233,*) "!the hrdat for 13.2 degs stack graphene"
   write(233,*) cot
   write(233,*) "9"
   write(233,*) "1   1   1   1   1   1   1   1   1"

   do m=-1,1
      do n=-1,1
         do i=1,cot
            do j=1,cot
				j1=atoms(i,:)
				j2=atoms(j,:)+n*va2+m*va1
               tm=hoppingvalue(j1,j2)
               tmc=hoppingvalue(atoms(i,:),atoms(j,:)+n*va2+m*va1)-hoppingvalue(atoms(j,:)+n*va2+m*va1,atoms(i,:))

               if(abs(tmc)>0.000001) stop
!~                if(i==1 .and. j==8) write(*,*) tm
               if(i==j .and. m==0 .and. n==0) then
                  write(233,*) m,n,0,i,j,'    ',0.0000,0.00000
               else
                  write(233,*) m,n,0,i,j,'    ',tm,0.00000
               end if
            end do
         end do
      end do
   end do
   write(*,*) 'hrokay'
end program

module ttms
   integer :: pjv(2)
   integer,parameter :: nlat=100000
contains
   function trans(v0,mat) result(v1)
      implicit none
      real :: v0(2),v1(2),mat(2,2),mativ(2,2)
      mativ=mat
      call inv_r(2,mativ)
      v1=matmul(mativ,v0)
      return
   end function
   function proj(m,theta) result(b)
      implicit none
      real :: theta,b(2),t1,t2,t3
      real :: m(2)
      real :: rang=3.1415926/3.0
      b=0
!~ 	m=m*316.0/246.0
      b(1)=m(1)*sin(rang-theta)/sin(rang)-m(2)*sin(theta)/sin(rang)
      b(2)=m(1)*sin(theta)/sin(rang)+m(2)*sin(2*rang-theta)/sin(rang)
      return
   end function proj
   subroutine kspc(n1,n2,kspace)
      implicit none
      integer :: n1,n2,i,j
      real :: kspace(n1*n2,2),tmmps
      do i=1,n1
         do j=1,n2
            kspace((i-1)*n1+j,1)=1.0*(2*i-n1-1)/(2*n1)
            kspace((i-1)*n1+j,2)=1.0*(2*j-n2-1)/(2*n2)
         end do
      end do
   end subroutine
   function dist(a,b) result(d)
      implicit none
      real :: a(2),b(2),d
      d=abs(a(1)-b(1))+abs(a(2)-b(2))
   end function

   function mindist(b) result(d)
      implicit none
      integer :: i
      real :: b(2),d,dn(7)
      real :: a(7,2)
      a(1,1)=floor(b(1))
      a(1,2)=floor(b(2))

      a(2,:)=a(1,:)+[1,0]
      a(3,:)=a(1,:)+[0,1]
      a(4,:)=a(1,:)+[1,1]
      a(5,:)=a(1,:)+[0.666,0.666]
      a(6,:)=a(5,:)+[-1.0,0.0]
      a(7,:)=a(5,:)+[0.0,-1.0]
      do i=1,7
         dn(i)=dist(a(i,:),b)
      end do
      d=minval(dn)
      do i=1,7
         if(dn(i)==d .and. d<0.01) then
            write(277,*) a(i,1),',',a(i,2)
            pjv=int(a(i,:))
         end if
      end do
!~ 	d(2)=0.2*(dn(1)+dn(2)+dn(3)+dn(4)+dn(5))
      return
   end function

   function gmat(length,lattice) result(mat) !find lattice vector from lattice
      implicit none
      real :: mat(2,2)
      integer :: length,lattice(length,2),absl(length),mxv
      integer,allocatable :: m(:),flgs(:)
      integer :: i,j,k
      
      do i=1,length
         absl(i)=abs(lattice(i,1))+abs(lattice(i,2))
      end do
      mxv=maxval(absl)
!~       write(*,*) mxv
      allocate(m(0:mxv))
      allocate(flgs(0:mxv))
      m=0
      do i=1,length
         m(absl(i))=i
      end do
      j=0
      do i=1,mxv
!~          if(j==1 .and. (mat(1,1)*lattice(m(i),2)-mat(2,1)*lattice(m(i),1))==0) cycle
         if(j>=2) go to 101
         if(m(i)/=0) then
            j=j+1
            mat(:,j)=lattice(m(i),:)
         end if
      end do
101   continue
   end function

end module



program main
   use ttms
   implicit none
   integer,parameter :: n1=255,n2=n1,nm=200,nn=600
   complex,parameter :: ii=(0.0,1.0)
   integer :: i,j,k,kn,lattice(nlat,2)
   real :: kspace(n1*n2,2),b(2),theta,b1(2),thres1
   real ::overlap,thres,m(2),ss,dsb,m1(2)
   real :: sb(nn),sbposition(nn,2)
   real :: mat(2,2),pjvec(2)

   real,parameter :: g1(2)=[2.46,0.00],g2(2)=[1.23,2.13]

   integer :: atoms,ai
   real,allocatable :: atompos(:,:)
   atoms=2
   allocate(atompos(atoms,2))
   atompos(1,:)=[0.333,0.333]
   atompos(2,:)=[0.666,0.666]

!~    call kspc(n1,n2,kspace)
!~ write(*,*) sin(1.57)

   open(UNIT=777,file='overlaprat.csv')
   open(unit=233,file='midpos.csv')
   open(unit=255,file='poscar.vasp')
   open(unit=277,file='layer1.csv')
   open(unit=234,file='layers.csv')


   kn=218

!~    read(*,*) mat(1,1),mat(2,1)
!~    read(*,*) mat(1,2),mat(2,2)

   sb=0
   sbposition=0
!~    do kn=1,nn

   overlap=0
   theta=real(kn)*3.1415926/1800.0
!~ write(*,*) theta,overlap
   ss=1.0/real(nm*nm)
!~ write(*,*) ss
!~    write(234,*) 0,',',0,',',0

   k=0
   lattice=0
   do i=1,nm
      do j=1,nm

         do ai=1,atoms
            m(1)=i-nm/2!real(i)
            m(2)=j-nm/2!real(j)
            m1=m+atompos(ai,:)-atompos(1,:)![0.666,0.666]
            b1=proj(m,theta)
            if(mindist(b1)<0.01  ) then
               k=k+1
               lattice(k,:)=pjv
               overlap=overlap+ss
               write(233,*) b1(1),',',b1(2)!i,',',j
            else
               overlap=overlap
            end if
         end do



         sb(kn)=sb(kn)+mindist(b)*ss
      end do
   end do
	write(*,*) k
   mat=gmat(nlat,lattice)
!~    mat=mat*4
!~    mat(:,1)=[1,2]*4
!~    mat(:,2)=[2,-3]*4
   write(*,*) mat(:,1)
   write(*,*) mat(:,2)
   write(234,*) '             '
   write(234,*) mat(:,1)
   write(234,*) mat(:,2)
   write(255,*) 'stack layers'
   write(255,*) 1.000000
   write(255,*) g1*mat(1,1)+mat(2,1)*g2,0.00
   write(255,*) g1*mat(1,2)+g2*mat(2,2),0.00
   write(255,*) 0.00,0.00,6.70
   write(255,*) 'C'
   write(255,*) 'Direct'
   k=0
   do i=1,nm
      do j=1,nm

         write(233,*) b(1),',',b(2)!i,',',j

         do ai=1,atoms
            m(1)=i-nm/2!real(i)
            m(2)=j-nm/2!real(j)
            m1=m+atompos(ai,:)![0.666,0.666]
            pjvec=trans(m1,mat)
            if(pjvec(1) > -0.01 .and. pjvec(1)<0.98 .and. pjvec(2) >-0.01 .and. pjvec(2) < 0.98) then
               write(234,*) 0,',',m1(1),',',m1(2)
               write(255,*) pjvec(1),pjvec(2),0.1
               k=k+1
            end if
            b1=proj(m1,theta)
            pjvec=trans(b1,mat)
            if(pjvec(1) > -0.01 .and. pjvec(1)<0.98 .and. pjvec(2) >-0.01 .and. pjvec(2) < 0.98) then
               write(234,*) 1,',',b1(1),',',b1(2)
               write(255,*) pjvec(1),pjvec(2),0.6
               k=k+1
            end if
            if(mindist(b1)<0.01  ) then
               overlap=overlap+ss
!~             write(233,*) b1(1),',',b1(2)!i,',',j

            else
               overlap=overlap
            end if
         end do

         sb(kn)=sb(kn)+mindist(b)*ss
      end do
   end do
!~    rewind(234)
   write(*,*) k
!~    write(234,*) k
   write(777,*) kn,',',overlap
!~    end do



end program main


 !> get the inverse of a real matrix
subroutine inv_r(ndim,Amat)

   implicit none

   integer,parameter :: dp=4

   integer           :: i
   integer           :: info

   integer,intent(in):: ndim

!    IPIV   : INTEGER. Array, DIMENSION at least max(1, n). The pivot indices that define
!    the permutation matrix P; row i of the matrix was interchanged with
!    row ipiv(i). Corresponds to the single precision factorization (if info=
!    0 and iter ≥ 0) or the double precision factorization (if info= 0 and
!    iter < 0).
   integer,allocatable   :: ipiv(:)


!    Amat  :
!    Overwritten by the factors L and U from the factorization of A = P*L*U;
!    the unit diagonal elements of L are not stored.
!    If iterative refinement has been successfully used (info= 0 and iter≥
!    0), then A is unchanged.
!    If double precision factorization has been used (info= 0 and iter <
!    0), then the array A contains the factors L and U from the factorization
!    A = P*L*U; the unit diagonal elements of L are not stored.
   real(dp),intent(inout):: Amat(ndim,ndim)

!    Bmat  :
!    Overwritten by the solution matrix X for dgesv, sgesv,zgesv,zgesv.
   real(dp),allocatable :: Bmat(:,:)


   allocate(ipiv(ndim))
   allocate(Bmat(ndim,ndim))

   ipiv=0

   ! unit matrix
   Bmat= 0d0
   do i=1,ndim
      Bmat(i,i)= 1d0
   enddo

   call sgesv(ndim,ndim,Amat,ndim,ipiv,Bmat,ndim,info)

   if(info.ne.0)print *,'something wrong with dgesv'

   Amat=Bmat

   return
end subroutine inv_r


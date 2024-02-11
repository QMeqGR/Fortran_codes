! Version
! 1.4.0  Wed Apr  8 10:16:34 CDT 2015
!       -add print_matrix funciont
!
! 1.3.0  Fri Feb 13 16:59:37 CST 2015
!       -add frac to cart function
!
! 1.2.2  Sun 18 Jan 2015
!       -add vtriple function
!
! NOTE:  Thu Jul 12 15:27:04 PDT 2007
!       -gfortran 4.1.1 produces code which gives
!        errors in the pbc routines, don't know why.
!        4.3.0 for sure does not.

! 1.2.1  Wed Oct  1 13:12:10 CDT 2008
!       -Fixed bug in matrix inverter for pivoting.
! 1.2.0  Tue Mar 20 14:40:17 PDT 2007
!       -Fixed bug in matrix inverter! This was a c-translation error
!        and the c routine works fine.
! 1.1.0  Mon Jan 22 20:01:13 PST 2007
! 1.0.0  Fri Oct 13 16:36:49 PDT 2006


!---------------------------------------
! Tolerance function
! Input:
!       a - input real
!       b - input real
!       tolerance - input real
!
! Output:
!       L - logical int
!
!---------------------------------------
subroutine tol(a,b,tolerance,same)

implicit none

logical, intent(out)             :: same
real*8, intent(in)               :: tolerance
real*8, intent(in)               :: a,b

if ( abs(b-a) .LT. tolerance ) then
   same=.true.
else
   same=.false.
end if

end subroutine tol

!---------------------------------------
! Distance between two points
!
! Input:
!       v1 - input vector
!       v2 - input vector
!
! Output:
!       d - distance
!
!---------------------------------------
subroutine dist(v1,v2,d)

implicit none

real*8, dimension(3)             :: v1,v2
real*8, intent(out)              :: d

d = sqrt( dot_product(v1-v2,v1-v2) )

end subroutine dist

!---------------------------------------
! Vector mag
!
! Input:
!       v - input vector
!
! Output:
!       d - vector magnitude
!
!---------------------------------------
subroutine vmag(v,d)

implicit none

real*8, dimension(3)             :: v
real*8, intent(out)              :: d

d = sqrt( dot_product(v,v) )

end subroutine vmag


!---------------------------------------
! Vector norm
!
! Input:
!       v - input vector
!
! Output:
!       w - norm of v
!
!---------------------------------------
subroutine vnorm(v,w)

implicit none

real*8, dimension(3)             :: v
real*8                           :: mag
real*8, dimension(3),intent(out) :: w

call vmag(v,mag)
w(:) = v(:)/mag

end subroutine vnorm


!--------------------------------------
! Cross product
! Input: a, b vectors
! Output: c vector
!--------------------------------------
subroutine cross(a,b,c)

  real*8, dimension(3), intent(in)     :: a,b
  real*8, dimension(3), intent(out)    :: c
  
  c(1) = a(2)*b(3)-a(3)*b(2)
  c(2) = -( a(1)*b(3)-a(3)*b(1) )
  c(3) = a(1)*b(2)-a(2)*b(1)

end subroutine cross


!---------------------------------------
! triple product for cell volume calc
! Input: input basis vectors
! Output: cell volume
!---------------------------------------
subroutine vtriple(ib,cvol)

real*8, dimension(3,3), intent(in) :: ib
real*8, dimension(3)               :: cprod
real*8, intent(out)                :: cvol

call cross(ib(1,:),ib(2,:),cprod)
cvol = dot_product(cprod,ib(3,:))

end subroutine vtriple

!---------------------------------------
! Transverse of a matrix function
! Input:  Mi
! Output: Mo
!---------------------------------------
subroutine trans33(Mi,Mo)

integer                             :: i
real*8, dimension(3,3), intent(in)  :: Mi
real*8, dimension(3,3), intent(out) :: Mo

do i=1,3
   Mo(:,i) = Mi(i,:)
end do

end subroutine trans33

!---------------------------------------
! Cartesian to Fractional coords
! Input: vi - input vector
!        ib - input basis (cart coords)
!     
!---------------------------------------
subroutine cart_to_frac(ib,vi,vo)

  implicit none

  real*8                            :: SMALL=1e-6
  real*8, dimension(3), intent(in)  :: vi
  real*8, dimension(3), intent(out) :: vo
  real*8, dimension(3,3)            :: ib,ib_inv,ib_inv_t

  ib_inv = ib
  call invrt_NxN(ib_inv,3)
  call trans33(ib_inv,ib_inv_t)
  call mat33_vec_mult( ib_inv_t , vi , vo )

  if ( vo(1) < -SMALL )   vo(1) = vo(1) + 1
  if ( vo(1) >  SMALL+1 ) vo(1) = vo(1) - 1
  if ( vo(2) < -SMALL )   vo(2) = vo(2) + 1
  if ( vo(2) >  SMALL+1 ) vo(2) = vo(2) - 1
  if ( vo(3) < -SMALL )   vo(3) = vo(3) + 1
  if ( vo(3) >  SMALL+1 ) vo(3) = vo(3) - 1

end subroutine cart_to_frac

!---------------------------------------
! Fractional to Cartesian coords
! Input: fi - input vector
!        ib - input basis (cart coords)
!     
!---------------------------------------
subroutine frac_to_cart(ib,fi,vo)

  implicit none

  real*8                            :: SMALL=1e-6
  real*8, dimension(3), intent(in)  :: fi
  real*8, dimension(3), intent(out) :: vo
  real*8, dimension(3,3)            :: ib,ib_t

  call trans33( ib , ib_t )
  call mat33_vec_mult( ib_t , fi , vo )

end subroutine frac_to_cart


!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!
! I have two routines for calculating distances with periodic
! boundary conditions.  The second one alters the location of
! the atoms!
!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!


!---------------------------------------
! Periodic boundary conditions
! Distance between two points
!
! Input:
!       v1 - input vector (cartesian coordinates)
!       v2 - input vector (cartesian coordinates)
!       ib - input basis (cartesian vectors)
!       ib_inv - input basis inverse
!
! Temp:
!       a,b,c - mag of lat vecs
!
! Output:
!       d - distance
!
!---------------------------------------
subroutine dist_pbc(v1,v2,ib,ib_inv,d)

implicit none

real*8, dimension(3)             :: v1,v2,dvec,dvecf
real*8, intent(out)              :: d
real*8, dimension(3,3)           :: ib,ib_inv,ib_inv_t

! decompose the vector between the two points of
! interest into components along the lattice vectors
dvec = v2-v1
call trans33(ib_inv,ib_inv_t)
call mat33_vec_mult( ib_inv_t , dvec , dvecf )

!print '(1X,3F10.3)', dvecf

if ( dvecf(1) >  0.5 ) dvec = dvec - ib(1,:)
if ( dvecf(1) < -0.5 ) dvec = dvec + ib(1,:)

if ( dvecf(2) >  0.5 ) dvec = dvec - ib(2,:)
if ( dvecf(2) < -0.5 ) dvec = dvec + ib(2,:)

if ( dvecf(3) >  0.5 ) dvec = dvec - ib(3,:)
if ( dvecf(3) < -0.5 ) dvec = dvec + ib(3,:)

d = sqrt(dot_product(dvec,dvec))

end subroutine dist_pbc

!---------------------------------------
! Periodic boundary conditions
! Distance between two points
! !!!!!!!!!!!!!!!!!!!!!!!!!!!!! ALTERS THE INPUT POSITIONS
!
! Input:
!       v1 - input vector (cartesian coordinates)
!       ib - input basis (cartesian vectors)
!       ib_inv - input basis inverse
!
! Input/Output:
!       v2 - input vector (cartesian coordinates)
! Temp:
!       a,b,c - mag of lat vecs
!
! Output:
!       d - distance
!
!---------------------------------------
subroutine dist_pbc_ALT(v1,v2,ib,d)

implicit none

integer                             :: debug=0
real*8, dimension(3)                :: v1,dvec,dvecf
real*8, dimension(3), intent(inout) :: v2
real*8, intent(out)                 :: d
real*8, dimension(3,3)              :: ib,ib_t_inv,ib_t

! decompose the vector between the two points of
! interest into components along the lattice vectors
! -for details see book: Project: Hydrides/MC, series T, book 4, p17
dvec = v2-v1
call trans33(ib,ib_t)
ib_t_inv = ib_t
call invrt_NxN(ib_t_inv,3)
call mat33_vec_mult( ib_t_inv , dvec , dvecf )

if ( debug==1 ) then
   print '(1X,1A10,3F10.3)', 'dvec =',dvec
   print '(1X,1A10,3F10.3)', 'dvecf=',dvecf
end if

if ( dvecf(1) >  0.5 ) then
   dvec = dvec - ib(1,:)
   v2 = v2 - ib(1,:)
end if
if ( dvecf(1) < -0.5 ) then
   dvec = dvec + ib(1,:)
   v2 = v2 + ib(1,:)
end if

if ( dvecf(2) >  0.5 ) then
   dvec = dvec - ib(2,:)
   v2 = v2 - ib(2,:)
end if
if ( dvecf(2) < -0.5 ) then
   dvec = dvec + ib(2,:)
   v2 = v2 + ib(2,:)
end if

if ( dvecf(3) >  0.5 ) then
   dvec = dvec - ib(3,:)
   v2 = v2 - ib(3,:)
end if
if ( dvecf(3) < -0.5 ) then
   dvec = dvec + ib(3,:)
   v2 = v2 + ib(3,:)
end if

if ( debug==1 ) then
   print '(1X,1A10,3F10.3)', 'ndvec =',dvec
   call mat33_vec_mult( ib_t_inv , dvec , dvecf )
   print '(1X,1A10,3F10.3)', 'ndvecf =',dvecf
end if

d = sqrt(dot_product(dvec,dvec))

end subroutine dist_pbc_ALT

!---------------------------------------
! Rotate by phi, theta, psi
!          (following Byron and Fuller)
! Input:
!       phi
!       the
!       psi
!       vi - input  vector
! Output:
!       vo - output vector
!---------------------------------------
subroutine rotate_ptp(phi,the,psi,vi,vo)
  
  implicit none
  
  real*8, dimension(3,3)             :: M
  real*8, dimension(3), intent(in)   :: vi
  real*8, dimension(3), intent(out)  :: vo
  real*8, intent(in)                 :: phi,the,psi
  real*8                             :: sphi,cphi,sthe,cthe,spsi,cpsi
  
  sphi=sin(phi); cphi=cos(phi)
  sthe=sin(the); cthe=cos(the)
  spsi=sin(psi); cpsi=cos(psi)
  
  M(1,1) =  cphi * cthe * cpsi - sphi * spsi;
  M(1,2) =  sphi * cthe * cpsi + cphi * spsi;
  M(1,3) = -sthe * cpsi;

  M(2,1) = -cphi * cthe * spsi - sphi * cpsi;
  M(2,2) = -sphi * cthe * spsi + cphi * cpsi;
  M(2,3) =  sthe * spsi;

  M(3,1) = cphi * sthe;
  M(3,2) = sphi * sthe;
  M(3,3) = cthe;

  call mat33_vec_mult(M,vi,vo)
  
end subroutine rotate_ptp


!---------------------------------------
! Multiply a vector by a 3x3 matrix
!
! Input:
!       Mt - matrix
!       vi - input  vector
!       vo - output vector
!
! Output: vo (new vector)
!
!---------------------------------------
subroutine mat33_vec_mult(Mt,vi,vo)

implicit none

integer                            :: i
real*8, dimension(3,3), intent(in) :: Mt
real*8, dimension(3), intent(in)   :: vi
real*8, dimension(3), intent(out)  :: vo



do i=1,3
   vo(i) = dot_product( Mt(i,:) , vi )
end do

end subroutine mat33_vec_mult



!---------------------------------------
! invert a simple NxN matrix
!
! Input:
!       Mt - matrix
!       N - rank of matrix
!
! Output: matrix Mt (inverted)
!
!---------------------------------------
subroutine invrt_NxN(Mt,n)

implicit none

logical                           :: same
integer                           :: debug=0
integer, intent(in)               :: n
integer                           :: i,j,l,k,r
real*8                            :: zero=0.0
real*8, dimension(n,n)            :: Mt
real*8, dimension(n,n)            :: Id
real*8, dimension(n)              :: tempm,tempi
real*8                            :: tolerance=1e-8


! create identity matrix
do i=1,n
   do j=1,n
      if ( i==j ) then
         Id(i,j) = 1.0
      else if ( i /= j ) then
         Id(i,j) = 0.0
      end if
   end do
end do

! the main loop
do l=1,n
   k=l
   ! in case diagonal element is 0
   call tol(Mt(k,l), zero, tolerance, same)
   if ( same ) then
      if ( debug > 0 ) print *, 'found diagonal 0 element'
      do r=1,n
         call tol(Mt(r,l), zero, tolerance, same)
         if ( same .eqv. .false. ) exit
      end do
      Mt(k,:) = Mt(k,:) + Mt(r,:)
      Id(k,:) = Id(k,:) + Id(r,:)
   end if

   ! get a 1 in row 1, col 1
   tempm(:) = Mt(k,:) / Mt(k,l)
   tempi(:) = Id(k,:) / Mt(k,l)
   Mt(k,:) = tempm(:)
   Id(k,:) = tempi(:)

   ! eliminate other elements in lth column
   do i=1,n
      if ( i==k .and. k==n ) exit
      if ( i==k .and. k/=n ) cycle
      ! mult kth row by ith column elem
      tempm(:) = Mt(k,:) * Mt(i,l)
      tempi(:) = Id(k,:) * Mt(i,l)
      ! perform the row op
      Mt(i,:) = Mt(i,:) - tempm(:)
      Id(i,:) = Id(i,:) - tempi(:)
   end do
end do

! replace Mt with Id
Mt=Id

end subroutine invrt_nxn

!====================================
!   print_matrix - Print matrix to stdout
!
!   Input:
!   M - matrix
!   nrows - number of rows
!   ncols - number of columns
!
!====================================
subroutine print_matrix(M,nrows,ncols,nfmt)
implicit none
character (len=32)              :: out_fmt
integer                         :: i, nfmt, nrows, ncols
real*8, allocatable, intent(inout) :: M(:,:)

! simple floats, no exponentials, increasing .xxx
if ( nfmt == 1 ) then
   write(out_fmt, '(a, i0, a)') '(',ncols,'F6.2)'
end if
if ( nfmt == 2 ) then
   write(out_fmt, '(a, i0, a)') '(',ncols,'F8.4)'
end if

! exponentials
if ( nfmt == 10 ) then
   write(out_fmt, '(a, i0, a)') '(',ncols,'E12.3e2)'
end if
if ( nfmt == 11 ) then
   write(out_fmt, '(a, i0, a)') '(',ncols,'E12.3e4)'
end if


if( allocated(M) ) then
   do i=1,nrows
      write (*,out_fmt) M(i,:)
   end do
end if

return
end subroutine print_matrix

!Program: anion_rotate
!
! This program takes two input files in the following formats:
!
! 1. Atom file format:
! (powder input file with lattice vectors
!  given in the file, as if -P option were given)
! (NAME MUST BE ar.dat)
!
! num_atoms
!
! c1x c1y c1z
! c2x c2y c2z  cartesian vecs (angstroms)
! c3x c3y c3z
!
! x1 y1 z1 Z1    FRACTIONAL atom coords and atomic number
! x2 y2 z2 Z2
! x3 y3 z3 Z3
! ...
! ...
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 2. Control file format: (NAME MUST BE control.dat)
!
! The format is listed and then the terms are defined:
!
! &condat
! debug = 1
! center_Z = 13
! vertex_Z =  1
! N_coord = 4
! dcv = 1.7
! derror = 0.3
! phi = 0.0
! the = 0.0
! psi = 0.0
! calc_symm = 1
! calc_dist = 1
! selective_rotate = 2
! an_rot(1) = 1
! an_rot(2) = 4
! /
!
! The terms are:
!
! debug_value          -- set to zero for 'powder' output
! center_Z vertex_Z N  -- atomic numbers of the anion center and vertices,
!                         anion coord (4 or 6)
! 1.23 0.05            -- anion-to-vertex length error ( dcv +/- error ) 
! phi theta psi        -- rotation angles defined as usual (Byron and Fuller)
!                         (input in DEGREES)
! n                    -- number of anions to selectively rotate
! n_1 n_2 ... n_m      -- list of anions to selectivley rotate
!                         (counting begins from 1)
! calc_symm            -- calculate the high symmetry directions
!                         for each of the anions and list them
!                         (1=on, 0=off)
! calc_dist            -- calculate the scale-invariant distortion
!                         for each of the anions and list them
!                         (1=on, 0=off)
!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! The output of the program will be identical to the input except that 
! the anions will be rotated by alpha,the,psia.  The output will
! be a 'powder' input file, usable with the -P switch.
!
!
! E. Majzoub
! Sandia National Labs
! 11 October 2006
!
! Version:
!
!  1.5.3  Wed Oct 24 10:24:14 CDT 2007
!        -fixed bug for 6-coord anions distortions.
!         The sigma part of the code needed the
!         same fix as 1.5.2 below for the d_ave part.
!  1.5.2  Wed Sep 12 14:27:34 CDT 2007
!        -fixed bug for 6-coord anions.
!         The distortion can't use all
!         the lengths since some are much
!         further apart.  See page 58
!         in MC book 4 (start july 2006).
!  1.5.1  Tue Mar 20 08:32:42 PDT 2007
!        -fixed bug in ehm_fortlib.f90
!  1.5.0  Mon Mar 19 09:05:26 PDT 2007
!        -uses namelist for input
!  1.4.0  Thu Mar 15 15:37:43 PDT 2007
!        -add calculation of distortion
!  1.3.0  Mon Jan 22 19:16:30 PST 2007
!        -add calculation of cryst directions
!         for each high symmetry axis of the
!         tetrahedron or octahedron
!  1.2.0  Mon Nov 13 14:01:27 PST 2006
!        -add support for selective
!         anion rotation
!  1.1.0  Fri Oct 13 16:45:35 PDT 2006
!        -add support for octahedrons
!  1.0.0  Fri Oct 13 16:36:49 PDT 2006
!        -works for tetrahedral anions


program anion_rotate

implicit none

! data structures
type anion
   integer                :: at_num
   real*8, dimension(7,3) :: vertex
end type anion

logical                            :: same,known_anion
integer                            :: nat_from_file
integer                            :: i,j,nat,ndat=0,inputstat,filein=10
integer                            :: k,l,m,sr,sr_flag
integer                            :: allocatestatus,debug=0
integer                            :: center_Z,vertex_Z,N_coord,num_anions=0
integer                            :: num_verts
integer                            :: selective_rotate=0
integer                            :: calc_symm=0
integer                            :: calc_dist=0, sigma_count=0
integer                            :: dist_count=0
integer                            :: condat_id=555
integer, dimension(100)            :: an_rot

real*8                             :: PI=3.14159265358979323846264
real*8                             :: R2D=57.2957795130823
real*8                             :: dcv,dij,derror,phi,the,psi
real*8                             :: d_ave,sigma
real*8                             :: t2mag,t3
real*8                             :: c2v_the,c2v_phi
real*8, dimension(3,3)             :: ib,ib_t,ib_t_inv,ib_inv,ib_temp
real*8, dimension(3)               :: t,t2
real*8, dimension(:,:),allocatable :: vir,vor,vic
integer, dimension(:),allocatable  :: Z,sel_rot
type(anion), dimension(:), allocatable :: an

namelist /condat/ debug, center_Z, vertex_Z, N_coord
namelist /condat/ dcv, derror, phi, the, psi
namelist /condat/ selective_rotate, calc_symm, calc_dist
namelist /condat/ an_rot

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!   READING input file ar.dat
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! count the number of lines in the file
open( filein, file='ar.dat',status='old',form='formatted',iostat=inputstat )
read (filein,*,iostat=inputstat) nat_from_file
do
   read (filein,*,iostat=inputstat) t(:)
!   read (*,*,iostat=inputstat) t(:)
   if (inputstat<0) exit
!   print '(I5,I5)',inputstat,ndat
   if ( inputstat==5008 ) exit !! 5008 is EOF
   ndat = ndat +1
end do

rewind filein
if ( debug>0 ) print '(1X,"ndat= ",i5)', ndat
nat = ndat-3
if ( debug>0 ) print '(1X,"nat = ",i5)', nat
if ( nat /= nat_from_file ) stop "**** error reading file ****"

! allocate memory for arrays
allocate(vir(ndat,3), stat=allocatestatus)
if ( allocatestatus /= 0 ) stop "*** not enough memory for vir ***"
allocate(vic(ndat,3), stat=allocatestatus)
if ( allocatestatus /= 0 ) stop "*** not enough memory for vic ***"
allocate(vor(ndat,3), stat=allocatestatus)
if ( allocatestatus /= 0 ) stop "*** not enough memory for vor ***"
allocate(Z(ndat), stat=allocatestatus)
if ( allocatestatus /= 0 ) stop "*** not enough memory Z ***"


! read in the cartesian vector input basis
read( filein,* ) nat_from_file
do i=1,3
   read( filein,* ) ib(i,:)
end do


! read in the atom positions
do i=1,nat
   read( filein,*, iostat=inputstat ) vir(i,:), Z(i)
!   read( *,*, iostat=inputstat ) vir(i,:)
!   print *,'inputstat=',inputstat
end do

close( filein )

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! READING CONTROL FILE control.dat
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
open(condat_id,file="control.dat",form="formatted")
read(condat_id, NML=condat)
close(condat_id)

allocate( sel_rot( selective_rotate ), stat=allocatestatus)
if ( allocatestatus /= 0 ) stop "*** not enough memory for sel_rot ***"

if ( debug>0 ) then
   print *, 'input parameters:'
   print '(1X,1A15,1I15)', "debug", debug
   print '(1X,1A15,1I15)', "ndat", ndat
   print '(1X,1A15,1I15)', "center_Z", center_Z
   print '(1X,1A15,1I15)', "vertex_Z", vertex_Z
   print '(1X,1A15,1I15)', "N_coord", N_coord
   print '(1X,1A15,1F15.10)', "dcv", dcv
   print '(1X,1A15,1F15.10)', "derror", derror
   print '(1X,1A15,1F15.10)', "phi", phi
   print '(1X,1A15,1F15.10)', "theta", the
   print '(1X,1A15,1F15.10)', "psi", psi
   print '(1X,1A15,1I15)', "selective_rotate", selective_rotate
   print '(1X,1A15,1I15)', "calc_symm", calc_symm
   print '(1X,1A15,1I15)', "calc_dist", calc_dist
   do i=1,selective_rotate
      sel_rot(i) = an_rot(i)
      print '(1X,1A15,1I15)', "rot anion no.", sel_rot(i)
   end do
   print *
end if

! convert the angles to radians
phi = phi*PI/180.0
the = the*PI/180.0
psi = psi*PI/180.0

! set the number of vertices (for loops)
if ( N_coord==4 ) num_verts=5
if ( N_coord==6 ) num_verts=7

! print the input basis
if ( debug>0 ) then
   print *, 'input (conventional) basis:'
   do i=1,3
      print '(1X,3F20.10)', ib(i,:)
   end do
   print *
end if

ib_temp = ib
call invrt_NxN(ib_temp,3)
call invrt_NxN(ib_temp,3)
if ( debug>0 ) then
   print *, 'inverse of inverse of input (conventional) basis:'
   do i=1,3
      print '(1X,3F20.10)', ib_temp(i,:)
   end do
   print *
end if



! print the input atoms
if ( debug>0 ) then
   print *, 'input atoms (fractional coordinates):'
   print '(1X,3A20,2A5)', "x","y","z","Z","i"
   do i=1,nat
      print '(1X,3F20.10,2I5)', vir(i,:), Z(i), i
   end do
   print *
end if


do i=1,nat
   if ( Z(i)==center_Z ) num_anions = num_anions + 1
end do
if ( debug>0 ) print *, 'number of anions should be = ', num_anions

! allocate memory for the anions
allocate( an(num_anions), stat=allocatestatus)
if ( allocatestatus /= 0 ) stop "*** not enough memory for anions ***"

! initialize the anion atom numbers to -1
! and the vetices to 1e10
do i=1,num_anions
   an(i)%at_num = -1
   do k=1,num_verts
      an(i)%vertex(k,1)=1e10
      an(i)%vertex(k,2)=1e10
      an(i)%vertex(k,3)=1e10
   end do
end do

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! program starts here
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! convert the input atoms to cartesian coordinates
do i=1,nat
   vic(i,:) = vir(i,1) * ib(1,:) + vir(i,2) * ib(2,:) + vir(i,3) * ib(3,:)
end do

! print the atoms in cartesian coordinates
if ( debug>0 ) then
   print *, 'Atoms in Cartesian coords:'
   print '(1X,3A20,2A5)', "x","y","z","Z","i"
   do i=1,nat
      print '(1X,3F20.10,2I5)', vic(i,:), Z(i), i
   end do
   print *
end if

! invert the transpose of the conv basis matrix (needed for dist_pbc function)
call trans33(ib,ib_t)
if ( debug > 1 ) then
   print *, 'transpose of conventional basis vectors:'
   do i=1,3
      print '(1X,3F20.10)', ib_t(i,:)
   end do
   print *
end if

! invert the transpose
ib_t_inv = ib_t
call invrt_NxN(ib_t_inv,3)

if ( debug > 2 ) then
   print *, 'inverse of transpose of conventional basis vectors:'
   do i=1,3
      print '(1X,3F20.10)', ib_t_inv(i,:)
   end do
   print *
end if

if ( debug > 2 ) then
   print *, 'Converting atoms back to reduced coords:'

   print '(1X,3A20,2A5)', "x","y","z","Z","i"
   do i=1,nat
      call mat33_vec_mult( ib_t_inv , vic(i,:) , t )
      print '(1X,3F20.10,2I5)', t, Z(i), i
   end do
end if

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! This section rotates the anions and requires several steps
! 1. Find the anions
!    a. find all atoms with Z=Zvertex at distance d+/-error (in p.b.c.)
!       (this number should equal the anion coordination)
!    b. save the atom coordinates to the particular anion to which they belong
! 2. Translate the anion to the origin
!    a. fix atoms that were on opposite side of cell (p.b.c. effect)
! 3. Rotate the anion
! 4. Translate back to it's initial position

if ( debug>0 ) then
   print *
   print *,' ======= Searching for anions ========'
   print *
   print '(1X,1A10,2A8,1A20,2A10)',&
        "","i","j","dij","Z(i)","Z(j)"
end if

! run through the anion centers
do i=1,nat
   if ( Z(i) /= center_Z ) cycle

   ! check to see if the center of this pair is already
   ! one of the anions
   known_anion=.false.
   do k=1,num_anions
      if ( an(k)%at_num == i ) known_anion=.true.
   end do
   if ( known_anion .eqv. .false. ) then
      ! find the first non negative at_num and place this as the new center
      do k=1,num_anions
         if ( an(k)%at_num == -1 ) then
            an(k)%at_num = i
            an(k)%vertex(1,:)=vic(i,:)
            if ( debug > 1 ) print '(1X,1A15,3F15.8,1I5)',&
                 'assigned CENTER', vic(i,:), k
            exit
         end if
      end do
   end if

   ! run through the anion vertices
   do j=1,nat
      if ( Z(j) /= vertex_Z ) cycle

!!    call dist_pbc( vic(i,:) , vic(j,:) , ib , dij )
      call dist_pbc_ALT( vic(i,:) , vic(j,:) , ib , dij )
      call tol(dij,dcv,derror,same)
      if ( debug > 2 ) print '(1X,1A10,2I8,1F20.5,2I10)', &
           "**dist**", i, j, dij, Z(i), Z(j)

      if ( same ) then
         if ( debug==1 ) print '(1X,1A10,2I8,1F20.12,2I10)',&
              'dist ',i,j,dij,Z(i),Z(j)
         do m=2,num_verts
            if ( an(k)%vertex(m,1) == 1e10 ) then
               an(k)%vertex(m,:) = vic(j,:)
               if ( debug > 1) print '(1X,1A15,3F15.8,1I5)',&
                    'assigned VERTEX', vic(j,:), m
               exit
            end if
         end do
            
      end if

   end do
end do


if ( debug>0 ) then
   print *
   print *,' ======= END searching for anions ========'
   print *
end if


! print the anions
if ( debug > 0 ) then
   print *,'Anions found (center is listed first)'
   do i=1,num_anions
      print '(1X,1A10,1I10)', "anion",i
      do k=1,num_verts
         print '(1X,3F20.10)', an(i)%vertex(k,:)
      end do
      print *
   end do
   
   ! c2v distances
   do i=1,num_anions
      print '(1X,1A10,1I10)', "anion",i
      do k=2,num_verts
         call dist_pbc(an(i)%vertex(1,:),an(i)%vertex(k,:),ib,ib_inv,dij)
         print '(1X,1A10,1F20.10)',"c2v dist=",dij
      end do
      print *
   end do
end if


!!!!!!!!!!!!!!!!!!!!!!
! TANSLATE ANIONS TO ORIGIN
if ( debug>0 ) then
   print *,'Translate the anions to the origin:'
end if

do i=1,num_anions
   do k=2,num_verts
      an(i)%vertex(k,:) = an(i)%vertex(k,:) - an(i)%vertex(1,:)
   end do
end do

! print the anions
if ( debug > 0 ) then
   print *,'(centers not translated)'
   do i=1,num_anions
      print '(1X,1A10,1I10)', "anion",i
      do k=2,num_verts
         print '(1X,3F20.10)', an(i)%vertex(k,:)
      end do
      print *
   end do
   
end if


!!!!!!!!!!!!!!!!!!!!!!
! ROTATE ANIONS ABOUT ORIGIN
if ( debug>0 ) then
   print *,'Rotate the anions:'
   if ( selective_rotate>0 ) print *, '(selective rotation)'
end if

do i=1,num_anions

   sr_flag=0
   do sr=1,selective_rotate
      if ( i == sel_rot(sr) ) sr_flag=1
   end do
   if ( sr_flag==0 .and. selective_rotate>0 ) cycle
   if ( debug>1 ) print *,'debug level 2: rotating anion ',i

   do k=2,num_verts
      
      t = an(i)%vertex(k,:)
      call rotate_ptp(phi,the,psi, t , an(i)%vertex(k,:) )
   end do
end do

! print the anions
if ( debug > 0 ) then
   do i=1,num_anions
      print '(1X,1A10,1I10)', "anion",i
      do k=2,num_verts
         print '(1X,3F20.10)', an(i)%vertex(k,:)
      end do
      print *
   end do
   
end if


!!!!!!!!!!!!!!!!!!!!!!
! TANSLATE ANIONS BACK
if ( debug>0 ) then
   print *,'Translate the anions back where they came from:'
end if

do i=1,num_anions
   do k=2,num_verts
      an(i)%vertex(k,:) = an(i)%vertex(k,:) + an(i)%vertex(1,:)
   end do
end do


! print the anions
if ( debug > 0 ) then
   print *,'(center is listed first)'
   do i=1,num_anions
      print '(1X,1A10,1I10)', "anion",i
      do k=1,num_verts
         print '(1X,3F20.10)', an(i)%vertex(k,:)
      end do
      print *
   end do
   
   ! c2v distances
   do i=1,num_anions
      print '(1X,1A10,1I10)', "anion",i
      do k=2,num_verts
         call dist_pbc(an(i)%vertex(1,:),an(i)%vertex(k,:),ib,ib_inv,dij)
         print '(1X,1A10,1F20.10)',"c2v dist=",dij
      end do
      print *
   end do
end if


! convert the positions to reduced coords for output
print '(1X,1I5)',nat
print *
do i=1,3
   print '(1X,3F20.15)', ib(i,:)
end do
print *
do i=1,nat
   if ( Z(i) /= center_Z .and. Z(i) /= vertex_Z ) then
      call cart_to_frac(ib,vic(i,:),t)
      print '(1X,3F20.15,1I5)', t, Z(i)
   end if
end do
print *
do i=1,num_anions
   do k=1,num_verts
      call cart_to_frac(ib,an(i)%vertex(k,:),t)
      if ( k==1 ) print '(1X,3F20.15,1I5)', t, center_Z
      if ( k>1 )  print '(1X,3F20.15,1I5)', t, vertex_Z
   end do
   print *
end do

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! calculate the symmetry directions of the anions if
! the calc_symm variable is set to 1. Do this AFTER
! the atom coordinates are converted to cartesian.
! Algorithm:
!   1. loop over all anions
!   2.   loop over all vertices for each anion
!   3.     calculate the vector ( vertex - center )
!   4.     normalize this vector
!   5.     calculate theta and phi values
!   6.     check to see if the dot prod of this vector
!          is within some tolerance of any vectors in
!          a list of these vectors, if not, add it
!          (to be implemented later)
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
if ( calc_symm == 1 ) then
   
   print *,'CALC_SYMM=1  normalized center to vertex vectors'
   do i=1,num_anions
      print '(1X,1A10,1I5,6A10)', "anion",i,"x","y","z","mag","theta","phi"
      do k=2,num_verts
         t = an(i)%vertex(k,:) - an(i)%vertex(1,:)
         call vnorm(t,t2)
         call vmag(t2,t2mag)
         c2v_the = acos( t2(3) )
         c2v_phi = atan( t2(2) / t2(1) )
         print '(1X,1A10,1I5,3F10.5,1F10.5,2F10.2)', &
              "vertx",k,t2(:),t2mag,c2v_the*R2D,c2v_phi*R2D
      end do
      print *
   end do

end if

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! calculate the scale-invariant distortion of the anions if
! the calc_dist variable is set to 1. Do this AFTER
! the atom coordinates are converted to cartesian.
! Algorithm:
!   1. loop over all anions
!   2.   loop over all vertices for each anion
!   3.     calc d_av = 1/6 sum_i<j mag( r_i - r_j )
!   4.     calc sigma = 1/d_av * ( sum_i<j ( d_av - mag( r_i-r_j ) )^2 )^(1/2)
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
if ( calc_dist == 1 ) then
   
   print *,'CALC_DIST=1  calculating scale-invariant distortions'
   do i=1,num_anions
      print '(1X,1A10,1I5,6A10)', &
           "anion",i,"vert k","vert l","mag","d_ave","sigma"

      d_ave = 0.0
      sigma = 0.0
      t3=0
      dist_count = 0
      sigma_count = 0

      do k=2,num_verts
         do l=2,k
            if ( l == k ) cycle
            t2 = an(i)%vertex(k,:) - an(i)%vertex(l,:)
            call vmag(t2,t2mag)
            if ( N_coord==6 .AND. (t2mag > (sqrt(2.0)*dcv + derror)) ) cycle
            print '(1X,1I25,1I10,1F10.5,2A10)', k,l,t2mag,"-","-"
            d_ave = d_ave + t2mag
            dist_count = dist_count + 1
         end do
      end do

      if ( num_verts==5 .AND. dist_count /= 6 ) stop "* distortion error *"
      if ( num_verts==7 .AND. dist_count /= 12 ) stop "* distortion error *"


      if ( num_verts == 5 ) d_ave = d_ave / 6.0;
      if ( num_verts == 7 ) d_ave = d_ave / 12.0;

      do k=2,num_verts
         do l=2,k
            if ( l == k ) cycle
            t2 = an(i)%vertex(k,:) - an(i)%vertex(l,:)
            call vmag(t2,t2mag)
            if ( N_coord==6 .AND. (t2mag > (sqrt(2.0)*dcv + derror)) ) cycle
            t3 = t3 + (d_ave - t2mag)**2
!            print '(1X,1F35.5,1F15.5)', t3, t2mag
            sigma_count = sigma_count + 1
         end do
      end do

      sigma = sqrt( t3/sigma_count ) / d_ave

      print '(1X,1F55.5,1F10.5)', d_ave,sigma
      print *

   end do
   
end if


end program anion_rotate



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

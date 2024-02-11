!Program: find_interstitials
!
! This program takes two input files in the following formats:
!
! 1. Atom file format:
! (powder input file with lattice vectors
!  given in the file, as if -P option were given)
! (NAME MUST BE fi_in.dat)
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
! N_coord = 4       ! 4...MAXCOORD, 4=tet 5=pent 6=oct
! dstep = 0.3       ! Angstrom, make larger for faster runs
! half_cell = 1     ! restricts to frac coords 0.0 to 0.5
!                     use this if you want to search a unit cell without repeats
!                     from symmetry ops when generating a full cell, for example
! iv_sig_min = 0.0 ! increase to find fewer interstitials
! iv_sig_max = 0.5 ! reduce to find fewer interstitials
! show_coord = 0   ! shows the coordinating atoms that make the site
! /
!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! The output of the program will be identical to the input except that 
! the anions will be rotated by alpha,the,psia.  The output will
! be a 'powder' input file, usable with the -P switch.
!
!
! E. Majzoub
! University of Missouri - St. Louis
! 18 Jan 2015
!
! Version:
!
!  1.0.0  Sun Jan 18 2015
!        -started from anion_rotate.f90 code v 1.5.3
!

#define MAXCOORD 12
#define MAXINTERS 10000

module data_types

  ! data structures
  type anion
     integer                :: at_num
     real*8, dimension(7,3) :: vertex
  end type anion
  
  type interstit
     integer                       :: num=0
     integer, dimension(MAXCOORD)  :: index=0
     real*8, dimension(MAXCOORD,3) :: vertex=0
     real*8, dimension(3)          :: center=0
     integer, dimension(MAXCOORD)  :: Zvert=0
     real*8                        :: d_ave=0,sigma=0,iv_sig=0,iv_ave=0
  end type interstit

end module data_types

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

program find_interstitials

use data_types
implicit none

integer                            :: nat_from_file,show_coord
integer                            :: i,k,nat,ndat=0,inputstat,filein=10
integer                            :: allocatestatus,debug=0
integer                            :: N_coord,n_pts
integer                            :: half_cell, NsampX, NsampY, NsampZ
integer                            :: cellnx=0,cellny=0,cellnz=0
integer                            :: condat_id=555
integer, dimension(3)              :: Nsamp
integer                            :: ix,iy,iz, itab
integer, dimension(MAXCOORD)       :: nn_table=0
integer                            :: shift_tab
integer                            :: MAX_INTERSTIT=MAXINTERS
integer, dimension(:),allocatable  :: Z
integer, dimension(3)              :: version
integer                            :: scmm=1

!real*8                             :: PI=3.14159265358979323846264
real*8                             :: dstep,dij
real*8                             :: sigmax=0,sigmin=1e9
real*8                             :: davemax=0,davemin=1e9
real*8                             :: cvol,tmp,iv_sig_max,iv_sig_min
real*8                             :: ivsigmin=1e9,ivsigmax=0
real*8                             :: ivavemin=1e9,ivavemax=0
real*8                             :: xtry=0,ytry=0,ztry=0
real*8, dimension(MAXCOORD  )      :: di_table=1e9
real*8, dimension(3,3)             :: ib,ib_t,ib_t_inv,ib_temp
real*8, dimension(3)               :: t,interstitial,pos_r,cart_tmp,cart_inter
real*8, dimension(MAXCOORD,3)      :: cart=0
real*8, dimension(:,:),allocatable :: vir,vor,vic

type(interstit), allocatable         :: inter(:)

namelist /condat/ debug, N_coord
namelist /condat/ dstep, half_cell,iv_sig_max
namelist /condat/ iv_sig_min, show_coord, scmm

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! VERSION
version(1) = 1
version(2) = 0
version(3) = 0
print '(1A10,1I3,1A1,1I3,1A1,1I3)','VERSION ', version(1), '.',version(2),'.',version(3)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!   READING input file fi_in.dat
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! count the number of lines in the file
open( filein, file='fi_in.dat',status='old',form='formatted',iostat=inputstat )
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
allocate(vir(nat,3), stat=allocatestatus)
if ( allocatestatus /= 0 ) stop "*** not enough memory for vir ***"
allocate(vic(nat,3), stat=allocatestatus)
if ( allocatestatus /= 0 ) stop "*** not enough memory for vic ***"
allocate(vor(nat,3), stat=allocatestatus)
if ( allocatestatus /= 0 ) stop "*** not enough memory for vor ***"
allocate(Z(nat), stat=allocatestatus)
if ( allocatestatus /= 0 ) stop "*** not enough memory Z ***"
allocate(inter(MAX_INTERSTIT), stat=allocatestatus)
if ( allocatestatus /= 0 ) stop "*** not enough memory for interstitials ***"

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

print *
print '(1X,1A5,1A30)', "", "====== control.dat ======="
print *
print '(1X,1A5,1A15,1I15)', "I:", "debug", debug
print '(1X,1A5,1A15,1I15)', "I:", "N_coord", N_coord
print '(1X,1A5,1A15,1I15)', "I:", "show_coord", show_coord
print '(1X,1A5,1A15,1F15.5)', "I:", "dstep", dstep
print '(1X,1A5,1A15,1F15.5)', "I:", "iv_sig_min", iv_sig_min
print '(1X,1A5,1A15,1F15.5)', "I:", "iv_sig_max", iv_sig_max
print '(1X,1A5,1A15,1I15)', "I:", "half_cell", half_cell
print '(1X,1A5,1A15,1I15)', "I:", "scmm", scmm
print *

! print the input basis
print *, 'input (conventional) basis:'
print *
do i=1,3
   print '(1X,1A5,3F20.10)', "IB:" , ib(i,:)
end do
print *


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


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! program starts here
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! print the cell volume
call vtriple(ib,cvol)
if ( debug>0 ) then
   print *, 'Cell volume:', cvol
end if

! Number of sample points for interstitials
n_pts = ceiling(cvol/dstep**3)
if ( debug>0 ) then
   print *, 'Number of sample points in cell', n_pts
end if

! calculate the sample points in each direction
call vmag(ib(1,:),tmp)
Nsamp(1) = ceiling(tmp/dstep)
call vmag(ib(2,:),tmp)
Nsamp(2) = ceiling(tmp/dstep)
call vmag(ib(3,:),tmp)
Nsamp(3) = ceiling(tmp/dstep)

if ( half_cell == 0 ) then
   NsampX = Nsamp(1); NsampY = Nsamp(2); NsampZ = Nsamp(3);
else
   NsampX = ceiling(real(Nsamp(1)/2));
   NsampY = ceiling(real(Nsamp(2)/2));
   NsampZ = ceiling(real(Nsamp(3)/2));      
end if

if ( debug>0 ) then
   print *, 'Sample points: ', Nsamp
end if

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

! This section searches for interstitial sites

print *
print *,' ======= interstitial sites ========'
print *

! run through the sample points
do ix=0, NsampX
   xtry = real(ix)/Nsamp(1)
   do iy=0, NsampY
      ytry = real(iy)/Nsamp(2)
      do iz=0, NsampZ
         ztry = real(iz)/Nsamp(3)

         interstitial(1)=xtry
         interstitial(2)=ytry
         interstitial(3)=ztry

         if ( debug > 1 ) then
            print '(1X,1A15,3F10.3)', 'Interstitial xyz: ', interstitial
         endif
         nn_table=0
         di_table=1e9


         do cellnx = -scmm,scmm
            do cellny = -scmm,scmm
               do cellnz = -scmm,scmm
                  do i=1,nat
                     pos_r(1) = cellnx + vir(i,1)
                     pos_r(2) = cellny + vir(i,2)
                     pos_r(3) = cellnz + vir(i,3)
                     
                     cart_tmp = pos_r(1) * ib(1,:) + pos_r(2) * ib(2,:) + pos_r(3) * ib(3,:)
                     cart_inter = interstitial(1) * ib(1,:) + &
                          interstitial(2) * ib(2,:) + &
                          interstitial(3) * ib(3,:)
                     call vmag(cart_inter - cart_tmp, dij)

                     if ( debug > 3 ) then
                        print '(1X,1A20,3I5,1I6,1F8.4)', 'nx ny nz i dij= ', cellnx, cellny, cellnz, i, dij
                     end if
                    
                     do itab=1,N_coord
                        if ( dij <= di_table(itab) ) then

                           if ( debug > 2 ) then
                              print *, 'New table entry, i itab', i, itab
                              print '(1X,3F10.4)', cart_tmp
                           end if
                           do shift_tab=N_coord,itab+1, -1
                              cart(shift_tab,:) = cart(shift_tab-1,:)
                              nn_table(shift_tab) = nn_table(shift_tab-1)
                              di_table(shift_tab) = di_table(shift_tab-1)
                              if ( debug > 3 ) then
                                 print '(1X,1A10,10I6)', "shft nn ", nn_table(:)
                                 print '(1X,1A10,10F6.3)', "shft di ", di_table(:)
                              end if
                           end do
                           nn_table(itab) = i
                           di_table(itab) = dij
                           cart(itab,:) = cart_tmp
                           if ( debug > 1 ) then
                              print '(1X,1A10,10I6)', "nn_tab= ", nn_table(:)
                              print '(1X,1A10,10F6.3)', "di_tab= ", di_table(:)
                           end if
                           exit
                        end if
                     end do
                     
                  end do
               end do
            end do
         end do
         if ( debug > 1 ) then
            print '(1X,1A10,10I6)', "nn_table= ", nn_table(:)
            print '(1X,1A10,10F6.3)', "di_table= ", di_table(:)
            do k=1,MAXCOORD
               print '(1X,1A8,1I2,3F10.4)','cart',k,cart(k,:)
            end do
         end if

         ! calc distortion and geometric center
         call get_sigma(inter,MAX_INTERSTIT,ib,cart,nn_table,N_coord,Z,nat,iv_sig_min,iv_sig_max,debug,show_coord)
         
         
      end do
   end do
end do

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  Get high and low sigma values
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
do i=1,MAX_INTERSTIT
   if ( inter(i)%Zvert(1) /= 0 .AND. inter(i)%sigma < sigmin ) then
      sigmin = inter(i)%sigma
   end if
   if ( inter(i)%Zvert(1) /= 0 .AND. inter(i)%sigma > sigmax ) then
      sigmax = inter(i)%sigma
   end if
   if ( inter(i)%Zvert(1) /= 0 .AND. inter(i)%d_ave < davemin ) then
      davemin = inter(i)%d_ave
   end if
   if ( inter(i)%Zvert(1) /= 0 .AND. inter(i)%d_ave > davemax ) then
      davemax = inter(i)%d_ave
   end if
   if ( inter(i)%Zvert(1) /= 0 .AND. inter(i)%iv_sig < ivsigmin ) then
      ivsigmin = inter(i)%iv_sig
   end if
   if ( inter(i)%Zvert(1) /= 0 .AND. inter(i)%iv_sig > ivsigmax ) then
      ivsigmax = inter(i)%iv_sig
   end if
   if ( inter(i)%Zvert(1) /= 0 .AND. inter(i)%iv_ave < ivavemin ) then
      ivavemin = inter(i)%iv_ave
   end if
   if ( inter(i)%Zvert(1) /= 0 .AND. inter(i)%iv_ave > ivavemax ) then
      ivavemax = inter(i)%iv_ave
   end if
end do
print *
print '(1X,1A10,1A15,1A15)',"","min","max"
print '(1X,1A10,2F15.6)', "sigma", sigmin,sigmax
print '(1X,1A10,2F15.6)', "d_ave", davemin,davemax
print '(1X,1A10,2F15.6)', "iv_ave", ivavemin,ivavemax
print '(1X,1A10,2F15.6)', "iv_sig", ivsigmin,ivsigmax


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  Print out fractional coordinates
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
if ( debug > 0 ) then
   print *
   print *, 'Basis'
   print *
   do i=1,3
      print '(1X,3F20.15)', ib(i,:)
   end do
end if
print *
print *, ' ========== Fractional coordinates =========='
print *
do i=1,MAX_INTERSTIT
   if ( inter(i)%Zvert(1) /= 0 ) then
      call cart_to_frac(ib,inter(i)%center,t)
      print '(1X,1A5,1I6,3F15.8)', "frac", i, t
   end if
end do
print *

end program find_interstitials


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine get_sigma(inter,MAX_INTERSTIT,ib,cart,nn_table,N_coord,Z,nat,iv_sig_min,iv_sig_max,debug,show_coord)

use data_types
implicit none

integer                                  :: dist_count,k,l,skip
integer, dimension(MAXCOORD), intent(in) :: nn_table
integer, intent(in)                      :: N_coord,nat,show_coord
integer                                  :: n_edges=0
integer, intent(in)                      :: debug,MAX_INTERSTIT
integer, save                            :: inter_count=0
integer, intent(in)                      :: Z(nat)

real*8                                   :: d_ave=0,sigma=0,t3=0,t2mag,iv_ave=0,iv_sig=0
real*8, intent(in)                       :: cart(MAXCOORD,3),iv_sig_max,iv_sig_min
real*8, dimension(3)                     :: t2,center=(/0,0,0/),center_frac,t_frac
real*8, dimension(3:3), intent(in)       :: ib
type(interstit)                          :: inter(MAX_INTERSTIT)

if ( debug > 1 ) then
   print *,'--- in get_sigma subroutine ---'
end if

! vert   lines
! 1      N-1
! 2      N-2
! 3      N-3
! ...
! N-1    N-(N-1)
! N      N-N=0
! sum these up
!   (N-1)*N -(1+2+3+...+(N-1))  =  N**2 - (N*(N+1)/2)
n_edges = N_coord**2 - (N_coord+1)*N_coord/2

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! calculate center, d_ave
center = (/0,0,0/)
dist_count=0; d_ave=0;
do k=1,N_coord
   if ( debug>2) then
      print '(1X,1A2,1I2,1A3,3F10.6)','--',k,'--',cart(k,:)
   end if
   center = center + cart(k,:)
   do l=1,k
      if ( l == k ) cycle
      t2 = cart(k,:) - cart(l,:)
      call vmag(t2,t2mag)
      if ( debug > 3 ) then
         print '(1X,1I25,1I10,1F10.5,2A10)', k,l,t2mag,"-","-"
      end if
      d_ave = d_ave + t2mag
      dist_count = dist_count + 1
   end do
end do
d_ave = d_ave / n_edges
center = center / N_coord
if ( debug > 2 ) then
   print '(1X,1A6,3F10.6)','--c--',center
end if

!!!!!!!!!!!!!!!!!!!!!!!!!
! calculate iv_ave, iv_sig
iv_ave=0; iv_sig=0;
do k=1,N_coord
   t2 = center - cart(k,:)
   call vmag(t2,t2mag)
   iv_ave = iv_ave + t2mag
end do
iv_ave = iv_ave / N_coord

do k=1,N_coord
   t2 = center - cart(k,:)
   call vmag(t2,t2mag)
   iv_sig = iv_sig + (iv_ave - t2mag )**2
end do
iv_sig = sqrt(iv_sig)/iv_ave


!!!!!!!!!!!!!!!!!!!!!!!!!
! calculate sigma
sigma=0; t3=0
do k=1,N_coord
   do l=1,k
      if ( l == k ) cycle
      t2 = cart(k,:) - cart(l,:)
      call vmag(t2,t2mag)
      t3 = t3 + (d_ave - t2mag)**2
   end do
end do
sigma = sqrt( t3 ) / d_ave

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  FIRST interstitial found
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
if ( inter_count == 0 .AND. iv_sig_min <= iv_sig .AND. iv_sig < iv_sig_max ) then
   inter(1)%num = 1
   inter(1)%center(:) = center
   inter(1)%d_ave = d_ave
   inter(1)%iv_ave = iv_ave
   inter(1)%iv_sig = iv_sig
   inter(1)%sigma = sigma
   inter(1)%index=nn_table
   do k=1,N_coord
      inter(1)%vertex(k,:)=cart(k,:)
      inter(1)%Zvert(k)=Z(nn_table(k))
   end do
   inter_count = 1
   print '(1X,1A9,1A6,1A9,1A9,1A9,1A8,1A8,1A8)', ' ', 'i', 'x','y','z','iv_ave','iv_sig','sigma'
   print *
   print '(1X,1A9,1I6,3F9.4,3F8.4)', 'inter', inter_count, center, iv_ave, iv_sig, sigma
   if ( debug > 0 .OR. show_coord == 1 ) then
      print *
      do k=1,N_coord
            t2 = center - cart(k,:)
            call vmag(t2,t2mag)
            print '(1X,1A9,1I6,3F9.4,1F8.4,1I6)', 'from', nn_table(k), cart(k,:),t2mag, Z(nn_table(k))
         end do
      print *
   end if
end if

! check to see if this center is already in the interstitial database
skip = 0;
if ( inter_count > 0 ) then
   do k=1,inter_count
      call cart_to_frac(ib,center,center_frac)
      call cart_to_frac(ib,inter(k)%center,t_frac)
      call vmag( t_frac - center_frac , t2mag)
      if ( t2mag < 0.0001 ) then
!         print *, 'k t2mag', k, t2mag
         skip = 1
         exit
      end if
   end do
   if ( skip == 0 .AND. iv_sig_min <= iv_sig .AND. iv_sig < iv_sig_max ) then
      inter_count = inter_count + 1
      print '(1X,1A9,1I6,3F9.4,3F8.4)', 'inter', inter_count, center, iv_ave, iv_sig, sigma
      if ( debug > 0 .OR. show_coord == 1) then
         print *
         do k=1,N_coord
            t2 = center - cart(k,:)
            call vmag(t2,t2mag)
            print '(1X,1A9,1I6,3F9.4,1F8.4,1I6)', 'from', nn_table(k), cart(k,:),t2mag, Z(nn_table(k))
         end do
         print *
      end if
      inter(inter_count)%num = 1
      inter(inter_count)%center(:) = center
      inter(inter_count)%d_ave = d_ave
      inter(inter_count)%iv_ave = iv_ave
      inter(inter_count)%iv_sig = iv_sig
      inter(inter_count)%sigma = sigma
      inter(inter_count)%index=nn_table
      do k=1,N_coord
         inter(inter_count)%vertex(k,:)=cart(k,:)
         inter(inter_count)%Zvert(k)=Z(nn_table(k))
      end do
      
   end if
end if

end subroutine get_sigma

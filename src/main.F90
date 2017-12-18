program main
! *******************************************************
! Fortran 90 driver for Fractal Meanfield Scattering Code
!

use fractal_meanfield_mod
use system_mod

implicit none


real(f), allocatable :: wavelength(:)   !! lambda in microns
real(f), allocatable :: rfi(:)          !! imaginary index of refraction for monomer ! for core if CORESHELL OPTION
real(f), allocatable :: rfr(:)          !! real index of refraction for monomer      ! for core if CORESHELL OPTION
real(f), allocatable :: nmonomer(:)     !! number of monomers
real(f), allocatable :: alpha(:)        !! packing coefficient
real(f), allocatable :: fractaldim(:)   !! fractal dimension
real(f), allocatable :: rmonomer(:)     !! monomer size (um) ! includes core+shell if CORESHELL OPTION
real(f), allocatable :: xv(:)           !! 
real(f), allocatable :: angle(:)        !! angle
real(f), allocatable :: rcore(:)        !! core radius (um)                        ! CORESHELL OPTION
real(f), allocatable :: rfiSH(:)        !! imaginary index of refraction for shell ! CORESHELL OPTION
real(f), allocatable :: rfrSH(:)        !! real index of refraction for shell      ! CORESHELL OPTION
real(f), allocatable :: lqext(:)        !! extinction efficiency
real(f), allocatable :: lqsca(:)        !! scattering efficiency
real(f), allocatable :: lasym(:)        !! asymmetry parameter
integer  ::  rc         !! error handling  
integer :: ok           !! input file checking
integer :: N, i
character(LEN=75) :: header
logical :: do_miess


write(*,*) "Running fractaloptics.exe ..."

!
! open input file
!
open(UNIT=1,FILE="INPUT",FORM="FORMATTED",status="OLD",ACTION="READ")
read(UNIT=1,FMT="(I5,L7)"),  N, do_miess
read(UNIT=1,FMT="(A75)"), header

!allocate arrays
allocate(wavelength(N))
allocate(rfi(N))
allocate(rfr(N))
allocate(nmonomer(N))
allocate(alpha(N))
allocate(fractaldim(N))
allocate(rmonomer(N))
allocate(xv(N))
allocate(angle(N))
allocate(rcore(N))
allocate(rfiSH(N))
allocate(rfrSH(N))
allocate(lqext(N))
allocate(lqsca(N))
allocate(lasym(N))

if (do_miess) then
  ! read input file with coreshell data see file "INPUT_CORESHELL"
  write(*,*) "--- coreshell calculation"  
  do, i=1,N
    read(UNIT=1,FMT="(f10.4,         f10.5,  f10.5,  f12.1,       f5.2,     f5.2,          f10.5,       f5.1,  f5.1,     f10.5,    f10.5,     f10.5)"), &
                      wavelength(i), rfi(i), rfr(i), nmonomer(i), alpha(i), fractaldim(i), rmonomer(i), xv(i), angle(i), rcore(i), rfiSH(i), rfrSH(i)
  enddo
else
  ! read input file for homogeneous monomer data see file "INPUT_MONOMER"
  write(*,*) "--- homogeneous monomer"  
  do, i=1,N
    read(UNIT=1,FMT="(f10.4,         f10.5,  f10.5,  f12.1,       f5.2,     f5.2,          f10.5,       f5.1,  f5.1)"), &
                      wavelength(i), rfi(i), rfr(i), nmonomer(i), alpha(i), fractaldim(i), rmonomer(i), xv(i), angle(i)
  enddo
  ! zero out coreshell data
  rcore(i) = 0.0_f
  rfiSH(i) = 0.0_f
  rfrSH(i) = 0.0_f
end if

close(UNIT=1)


do i=1,N
  write(*,*) "   computing INPUT entry ",i
  call fractal_meanfield(do_miess, &          !! logical do coreshell calculation
                         wavelength(i), &     !! lambda in microns
                         rfi(i), &            !! imaginary index of refraction (core)
                         rfr(i), &            !! real index of refraction (core)
                         nmonomer(i), &       !! number of monomers
                         alpha(i), &          !! packing coefficient
                         fractaldim(i), &     !! fractal dimension 
                         rmonomer(i), &       !! monomer size (total = core+shell)
                         xv(i), &             !! xv,"set to 1"
                         angle(i), &          !! angle, set to 0
                         rcore(i), &          !! core radius in microns
                         rfiSH(i), &          !! imaginary index of refraction (shell)
                         rfrSH(i), &          !! real index of refraction (shell)
                         lqext(i), &          !! extinction efficiency
                         lqsca(i), &          !! scattering efficiency
                         lasym(i), &          !! asymmetry parameter
                         rc)

enddo

!
! open output file
!
open(UNIT=2,FILE="OUTPUT",FORM="FORMATTED",status="REPLACE",ACTION="WRITE")
write(UNIT=2,FMT="(I5)"),  N
write(UNIT=2,FMT="(5A15)"), "WVL(um)", "QEXT","QSCA","SINGSCAT","ASYM"
do i=1,N
  write(UNIT=2,FMT="(f15.7,         f15.7,    f15.7,    f15.7,             f15.7)"), &
                     wavelength(i), lqext(i), lqsca(i), lqsca(i)/lqext(i), lasym(i)  
enddo
close(UNIT=2)


!write(*,*) "Qext", lqext
!write(*,*) "Qsca", lqsca
!write(*,*) "SSA",  lqsca/lqext
!write(*,*) "Asym", lasym
write(*,*) "Data written successfully to OUTPUT"


!deallocate arrays
deallocate(wavelength)
deallocate(rfi)
deallocate(rfr)
deallocate(nmonomer)
deallocate(alpha)
deallocate(fractaldim)
deallocate(rmonomer)
deallocate(xv)
deallocate(angle)
deallocate(rcore)
deallocate(rfiSH)
deallocate(rfrSH)
deallocate(lqext)
deallocate(lqsca)
deallocate(lasym)

end program main
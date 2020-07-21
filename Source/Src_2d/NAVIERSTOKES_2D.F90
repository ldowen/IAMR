
#undef  BL_LANG_CC
#ifndef BL_LANG_FORT
#define BL_LANG_FORT
#endif
  
#include <AMReX_REAL.H>
#include <AMReX_CONSTANTS.H>
#include <AMReX_BC_TYPES.H>
#include <NAVIERSTOKES_F.H>
#include <AMReX_ArrayLim.H>

#define SDIM 2


module navierstokes_2d_module
  
  implicit none

  private 

  public gradp, fort_putdown, fort_maxval, &
       summass, summass_eb, cen2edg
  
contains

      subroutine gradp (&
          p,DIMS(p),&
          gp,DIMS(gp),&
          lo,hi,dx,is_full) bind(C,name="gradp")
!c 
!c     Compute a cell centered gradient from a node
!c     centered field.  Returns all components of GRADP
!c
      implicit none
      integer    DIMDEC(p)
      integer    DIMDEC(gp)
      integer    lo(SDIM), hi(SDIM)
      REAL_T     dx(SDIM)
      REAL_T     p(DIMV(p))
      REAL_T     gp(DIMV(gp),SDIM)
      integer    is_full
      integer    i,j
      REAL_T     ddx, ddy

      write(*,*) " in gradp"

      if (is_full .eq. 0) then
         ddx = half/dx(1)
         ddy = half/dx(2)

        do j = lo(2), hi(2)
        do i = lo(1), hi(1)

          gp(i,j,1) = ddx*(p(i+1,j)-p(i,j)+p(i+1,j+1)-p(i,j+1))
          gp(i,j,2) = ddy*(p(i,j+1)-p(i,j)+p(i+1,j+1)-p(i+1,j))
        
   ! Below is the same "stencil" as the compute fluxes in AMReX     
!          gp(i,j,1) = ddx*(-p(i,j)+p(i+1,j)-p(i,j+1)+p(i+1,j+1))
!          gp(i,j,2) = ddy*(-p(i,j)-p(i+1,j)+p(i,j+1)+p(i+1,j+1))

        end do
        end do
      else
        do j = lo(2), hi(2)
        do i = lo(1), hi(1)
          gp(i,j,1) = (p(i+1,j)-p(i,j)+p(i+1,j+1)-p(i,j+1))
          gp(i,j,2) = (p(i,j+1)-p(i,j)+p(i+1,j+1)-p(i+1,j))
        end do
        end do
      endif

    end subroutine gradp

!c :: ----------------------------------------------------------
!c :: Replace coarse grid pressure data with corresponding
!c :: fine grid pressure data.
!c ::
!c :: INPUTS / OUTPUTS:
!c ::  crse      <=  coarse grid data
!c ::  DIMS(crse) => index limits of crse
!c ::  fine       => fine grid data
!c ::  DIMS(fine) => index limits of fine 
!c ::  lo,hi      => index limits of overlap (crse grid)
!c ::  ratios     => IntVect refinement ratio
!c ::
!c :: NOTE:
!c ::  Assumes pressure fields are node based
!c :: ----------------------------------------------------------
!c ::
      subroutine fort_putdown (crse,DIMS(crse),&
           fine,DIMS(fine),lo,hi,ratios) &
           bind(C,name="fort_putdown")
      implicit none
      integer  DIMDEC(crse)
      integer  DIMDEC(fine)
      integer  lo(2), hi(2)
      integer  ratios(2)
      REAL_T   crse(DIMV(crse))
      REAL_T   fine(DIMV(fine))

      integer  ic, jc
      integer  lratx, lraty

      lratx = ratios(1)
      lraty = ratios(2)

      do jc = lo(2), hi(2)
         do ic = lo(1), hi(1)
            crse(ic,jc) = fine(lratx*ic,lraty*jc)
         end do
      end do

      write(*,*) " in fort_putdown"
    end subroutine fort_putdown

!c :: ----------------------------------------------------------
!c :: MAXVAL
!c ::             maxval = max{ rho(i,j) }
!c ::
!c :: ----------------------------------------------------------
!c ::
     subroutine fort_maxval(rho,DIMS(rho),DIMS(grid),mxval) &
          bind(C,name="fort_maxval")

       implicit none
       integer DIMDEC(rho)
       integer DIMDEC(grid)
       REAL_T  rho(DIMV(rho))
       REAL_T  mxval

       integer i,j

       mxval = -Huge(zero)

       do j = ARG_L2(grid), ARG_H2(grid)
          do i = ARG_L1(grid), ARG_H1(grid)
             mxval = max(mxval, rho(i,j))
          end do
       end do
      write(*,*) " in fort_maxval"

     end subroutine fort_maxval

!c :: ----------------------------------------------------------
!c :: SUMMASS
!c ::             MASS = sum{ vol(i,j)*rho(i,j) }
!c ::
!c :: INPUTS / OUTPUTS:
!c ::  rho        => density field
!c ::  DIMS(rho)  => index limits of rho aray
!c ::  lo,hi      => index limits of grid interior
!c ::  dx	 => cell size
!c ::  mass      <=  total mass
!c ::  r		 => radius at cell center
!c ::  irlo,hi    => index limits of r array
!c ::  rz_flag    => == 1 if R_Z coords
!c ::  tmp        => temp column array
!c :: ----------------------------------------------------------
!c ::
       subroutine summass(rho,DIMS(rho),DIMS(grid),dx,mass,&
            r,irlo,irhi,rz_flag) bind(C,name="summass")

       implicit none
       integer irlo, irhi, rz_flag
       integer DIMDEC(rho)
       integer DIMDEC(grid)
       REAL_T  mass, dx(2)
       REAL_T  rho(DIMV(rho))
       REAL_T  r(irlo:irhi)

       integer i, j
       REAL_T  dr, dz, vol

       dr = dx(1)
       dz = dx(2)

       mass = zero

       do i = ARG_L1(grid), ARG_H1(grid)
          vol = dr*dz
	  if (rz_flag .eq. 1) vol = vol*two*Pi*r(i)
          do j = ARG_L2(grid), ARG_H2(grid)
	     mass = mass + vol*rho(i,j)
	  end do
       end do
      write(*,*) " in summass"

     end subroutine summass

!c :: ----------------------------------------------------------
!c :: SUMMASS
!c ::             MASS = sum{ volfrac*vol(i,j)*rho(i,j) }
!c ::
!c :: INPUTS / OUTPUTS:
!c ::  rho        => density field
!c ::  DIMS(rho)  => index limits of rho aray
!c ::  lo,hi      => index limits of grid interior
!c ::  dx	 => cell size
!c ::  mass      <=  total mass
!c ::  r		 => radius at cell center
!c ::  irlo,hi    => index limits of r array
!c ::  rz_flag    => == 1 if R_Z coords
!c ::  tmp        => temp column array
!c :: ----------------------------------------------------------
!c ::
     subroutine summass_eb(rho,DIMS(rho),DIMS(grid),vf,DIMS(vf),dx,mass,&
            r,irlo,irhi,rz_flag) bind(C,name="summass_eb")

       implicit none
       integer irlo, irhi, rz_flag
       integer DIMDEC(rho)
       integer DIMDEC(vf)
       integer DIMDEC(grid)
       REAL_T  mass, dx(2)
       REAL_T  rho(DIMV(rho))
       REAL_T  vf(DIMV(vf))
       REAL_T  r(irlo:irhi)

       integer i, j
       REAL_T  dr, dz, vol

       dr = dx(1)
       dz = dx(2)

       mass = zero

       do i = ARG_L1(grid), ARG_H1(grid)
          vol = dr*dz
	  if (rz_flag .eq. 1) vol = vol*two*Pi*r(i)
          do j = ARG_L2(grid), ARG_H2(grid)
	     mass = mass + vf(i,j)*vol*rho(i,j)
	  end do
       end do
      write(*,*) " in summass_eb"

     end subroutine summass_eb
     
!c ::
!c :: ----------------------------------------------------------
!c :: This routine fills an edge-centered fab from a cell-centered
!c :: fab using simple linear interpolation.
!c ::
!c :: INPUTS / OUTPUTS:
!c ::  lo,hi      => index limits of the of the cell-centered fab
!c ::  DIMS(cfab) => index limits of the cell-centered fab
!c ::  cfab       => cell-centered data
!c ::  DIMS(efab) => index limits of the edge-centered fab
!c ::  efab       => edge-centered fab to fill
!c ::  n!c         => Number of components in the fab to fill
!c ::  dir        => direction data needs to be shifted to get to edges
!c :: ----------------------------------------------------------
!c ::
      subroutine cen2edg(lo, hi, &
          DIMS(cfab), cfab,&
          DIMS(efab), efab, nc, dir,&
          isharm) bind(C,name="cen2edg")
      implicit none
      integer lo(SDIM), hi(SDIM), nc, dir, isharm
      integer DIMDEC(cfab)
      integer DIMDEC(efab)
      REAL_T  cfab(DIMV(cfab), nc)
      REAL_T  efab(DIMV(efab), nc)

      integer i,j,n

      if ( isharm .eq. 0 ) then
         if (dir .EQ. 0) then
            do n = 1,nc
               do j = lo(2), hi(2)
                  do i = lo(1), hi(1)
                     efab(i,j,n) = half*(cfab(i,j,n) + cfab(i-1,j,n))
                  end do
               end do
            end do
         else
            do n = 1,nc
               do j = lo(2), hi(2)
                  do i = lo(1), hi(1)
                     efab(i,j,n) = half*(cfab(i,j,n) + cfab(i,j-1,n))
                  end do
               end do
            end do
         end if
      else
         if (dir .EQ. 0) then
            do n = 1,nc
               do j = lo(2), hi(2)
                  do i = lo(1), hi(1)
                     if((cfab(i,j,n) * cfab(i-1,j,n)).gt.zero)then
                        efab(i,j,n)&
                            = 2*(cfab(i,j,n) * cfab(i-1,j,n))/&
                            (cfab(i,j,n) + cfab(i-1,j,n))
                     else
                        efab(i,j,n)=zero
                     endif
                  end do
               end do
            end do
         else
            do n = 1,nc
               do j = lo(2), hi(2)
                  do i = lo(1), hi(1)
                     if((cfab(i,j,n) * cfab(i,j-1,n)).gt.zero)then
                        efab(i,j,n)&
                            = 2*(cfab(i,j,n) * cfab(i,j-1,n))/&
                            (cfab(i,j,n) + cfab(i,j-1,n))
                     else
                        efab(i,j,n)=zero
                     endif
                  end do
               end do
            end do
         end if
      end if
      write(*,*) " in cen2edg"
    end subroutine cen2edg

  end module navierstokes_2d_module

module mrg_x2{c}_mct
!<list> key="use" number=1
  use shr_kind_mod
  use mct_mod
  use seq_flds_mod
  use seq_flds_indices
  use seq_comm_mct
  use seq_cdata_mod
!</list>

!<list> key="use" number=2
  use shr_kind_mod, only: r8 => shr_kind_r8, cl => shr_kind_cl
  use mct_mod
  use seq_flds_mod
  use seq_flds_indices
  use seq_comm_mct
  use seq_cdata_mod
  use seq_infodata_mod
!</list>

  implicit none
  save
  private ! except

!--------------------------------------------------------------------------
! Public interfaces
!--------------------------------------------------------------------------

  public :: mrg_x2{c}_init_mct
  public :: mrg_x2{c}_run_mct
  public :: mrg_x2{c}_final_mct

!--------------------------------------------------------------------------
! Private data
!--------------------------------------------------------------------------

!===========================================================================================
contains
!===========================================================================================

!<list> key="cother" number=1
  subroutine mrg_x2{c}_init_mct( cdata_{c} ,l2x_{c}, o2x_{c}, i2x_{c}, x{c}o_{c} ) 
!</list>

!<list> key="cother" number=3
  subroutine mrg_x2{c}_init_mct( cdata_{c} ,a2x_{c}, i2x_{c}, r2x_{c})
!</list>

!<list> key="cother" number=2
  subroutine mrg_x2{c}_init_mct( cdata_{c} ,a2x_{c} )
!</list>

!<list> key="cother" number=4
  subroutine mrg_x2{c}_init_mct( cdata_{c} ,a2x_{c}, o2x_{c} )
!</list>

!<list> key="cother" number=5
  subroutine mrg_x2{c}_init_mct( cdata_{c} ,s2x_{c} )
!</list>

!<list> key="cother" number=6
  subroutine mrg_x2{c}_init_mct( cdata_{c} ,g2x_{c} )
!</list>

    type(seq_cdata) ,intent(in)     :: cdata_{c}
!<list> key="cother" number=1
    type(mct_aVect), intent(inout) :: l2x_{c}
    type(mct_aVect), intent(inout) :: o2x_{c}
    type(mct_aVect), intent(inout) :: i2x_{c}
    type(mct_aVect), intent(inout) :: x{c}o_{c}
!</list>

!<list> key="cother" number=2
    type(mct_aVect), intent(inout)  :: a2x_{c}
!</list>

!<list> key="cother" number=3
    type(mct_aVect), intent(inout)  :: a2x_{c}
    type(mct_aVect), intent(inout)  :: i2x_{c}
    type(mct_aVect), intent(inout)  :: r2x_{c}
!</list>

!<list> key="cother" number=4
    type(mct_aVect), intent(inout)  :: a2x_{c}
    type(mct_aVect), intent(inout)  :: o2x_{c}
!</list>

!<list> key="cother" number=5
    type(mct_aVect), intent(inout)  :: s2x_{c}
!</list>

!<list> key="cother" number=6
    type(mct_aVect), intent(inout)  :: g2x_{c}
!</list>

    type(mct_GsMap), pointer        :: GSMap_{ccc}
    integer                         :: mpicom
!<list> key="def_var" number=1
    integer                  :: lsize
!</list>

    ! Set gsMap
    call seq_cdata_setptrs(cdata_{c}, gsMap=gsMap_{ccc}, mpicom=mpicom)

!<list> key="def_var" number=1
    lsize = mct_GSMap_lsize(GSMap_{ccc}, mpicom)
!</list>

!<list> key="cother" number=1
    ! Initialize av for land export state on atmosphere grid/ decomp

    call mct_aVect_init(l2x_{c}, rList=seq_flds_l2x_fields, lsize=lsize)
    call mct_aVect_zero(l2x_{c})

    ! Initialize av for ocn export state on atmosphere grid/decomp

    call mct_aVect_init(o2x_{c}, rList=seq_flds_o2x_fields, lsize=lsize)
    call mct_aVect_zero(o2x_{c})

    ! Initialize av for ice export state on atmosphere grid/decomp

    call mct_aVect_init(i2x_{c}, rList=seq_flds_i2x_fields, lsize=lsize)
    call mct_aVect_zero(i2x_{c})

    ! Initialize av for atm/ocn flux calculation on atmosphere grid/decomp

    call mct_aVect_init(x{c}o_{c}, rList=seq_flds_x{c}o_fields, lsize=lsize)
    call mct_aVect_zero(x{c}o_{c})
!</list>

!<list> key="cother" number=2
    ! Initialize av for atmosphere export state on lnd decomp

    call mct_aVect_init(a2x_{c}, rList=seq_flds_a2x_fields, &
         lsize=mct_gsMap_lsize(gsMap_{ccc}, mpicom))
    call mct_aVect_zero(a2x_{c})
!</list>

!<list> key="cother" number=3
    ! Initialize av for atmosphere export state on ocn decomp

    call mct_aVect_init(a2x_{c}, rList=seq_flds_a2x_fields, &
               lsize=MCT_GSMap_lsize(GSMap_{ccc}, mpicom))
    call mct_aVect_zero(a2x_{c})

    ! Initialize av for ice export state on ocn decomp

    call mct_aVect_init(i2x_{c}, rList=seq_flds_i2x_fields, &
               lsize=mct_GSMap_lsize(GSMap_{ccc}, mpicom))
    call mct_aVect_zero(i2x_{c})

    ! Initialize av for rof export state on ocn decomp

    call mct_aVect_init(r2x_{c}, rList=seq_flds_r2x_fields, &
               lsize=mct_GSMap_lsize(GSMap_{ccc}, mpicom))
    call mct_aVect_zero(r2x_{c})
!</list>

!<list> key="cother" number=4
    ! Initialize av for atmosphere export state on ice decomp

    call MCT_aVect_init(a2x_{c}, rList=seq_flds_a2x_fields, &
               lsize=mct_gsMap_lsize(gsMap_{ccc}, mpicom))
    call MCT_aVect_zero(a2x_{c})

    ! Initialize av for ocn export state on ice decomp

    call MCT_aVect_init(o2x_{c}, rList=seq_flds_o2x_fields, &
               lsize=mct_gsMap_lsize(gsMap_{ccc}, mpicom))
    call MCT_aVect_zero(o2x_{c})
!</list>

!<list> key="cother" number=6
    ! Initialize av for atmosphere export state on lnd decomp

    call mct_aVect_init(g2x_{c}, rList=seq_flds_g2x_fields, &
         lsize=mct_gsMap_lsize(gsMap_{ccc}, mpicom))
    call mct_aVect_zero(g2x_{c})
!</list>

!<list> key="cother" number=5
    ! Initialize av for atmosphere export state on lnd decomp

    call mct_aVect_init(s2x_{c}, rList=seq_flds_s2x_fields, &
         lsize=mct_gsMap_lsize(gsMap_{ccc}, mpicom))
    call mct_aVect_zero(s2x_{c})
!</list>

  end subroutine mrg_x2{c}_init_mct

!===========================================================================================
!<list> key="cother2" number=1
  subroutine mrg_x2{c}_run_mct( cdata_{c}, l2x_{c}, o2x_{c}, x{c}o_{c}, i2x_{c}, fractions_{c}, x2{c}_{c} ) 
!</list>

!<list> key="cother2" number=2
  subroutine mrg_x2{c}_run_mct( cdata_{c}, a2x_{c}, x2{c}_{c} )
!</list>

!<list> key="cother2" number=3
  subroutine mrg_x2{c}_run_mct( cdata_{c}, a2x_{c}, i2x_{c}, xa{c}_{c}, fractions_{c}, x2{c}_{c} )
!</list>

!<list> key="cother2" number=4
  subroutine mrg_x2{c}_run_mct( cdata_{c}, a2x_{c}, o2x_{c}, x2{c}_{c} )
!</list>

!<list> key="cother2" number=5
  subroutine mrg_x2{c}_run_mct( cdata_{c}, s2x_{c}, x2{c}_{c} )
!</list>

!<list> key="cother2" number=6
  subroutine mrg_x2{c}_run_mct( cdata_{c}, g2x_{c}, x2{c}_{c} )
!</list>

    !----------------------------------------------------------------------- 
    !
    ! Arguments
    !
    type(seq_cdata), intent(in)     :: cdata_{c}
!<list> key="cother2" number=1
    type(mct_aVect), intent(in)     :: l2x_{c}
    type(mct_aVect), intent(in)     :: o2x_{c}
    type(mct_aVect), intent(in)     :: x{c}o_{c}
    type(mct_aVect), intent(in)     :: i2x_{c}
    type(mct_aVect), intent(inout)  :: x2{c}_{c}
!</list>

!<list> key="cother2" number=2
    type(mct_aVect), intent(inout)     :: a2x_{c}
    type(mct_aVect), intent(inout)  :: x2{c}_{c}
!</list>

!<list> key="cother2" number=3
    type(mct_aVect), intent(in)     :: a2x_{c}
    type(mct_aVect), intent(in)     :: xa{c}_{c}
    type(mct_aVect), intent(in)     :: i2x_{c}
    type(mct_aVect), intent(inout)  :: x2{c}_{c}
!</list>

!<list> key="cother2" number=4
    type(mct_aVect), intent(in)     :: a2x_{c}
    type(mct_aVect), intent(in)     :: o2x_{c}
    type(mct_aVect), intent(inout)  :: x2{c}_{c}
!</list>

!<list> key="cother2" number=5
    type(mct_aVect), intent(in)     :: s2x_{c}
    type(mct_aVect), intent(inout)  :: x2{c}_{c}
!</list>

!<list> key="cother2" number=6
    type(mct_aVect), intent(in)     :: g2x_{c}
    type(mct_aVect), intent(inout)  :: x2{c}_{c}
!</list>

!<list> key="fractions" number=1
    type(mct_aVect), intent(in)     :: fractions_{c}
!</list>

    !
    ! Local variables
    !
    logical :: usevector    ! use vector-friendly mct_copy
!<list> key="run_def" number=1
    integer :: i
    real(r8):: flux_epbalfact
    character(len=cl) :: flux_epbal
    type(seq_infodata_type),pointer :: infodata
    !
    !----- formats -----
    character(*),parameter :: F01 =   "('(mrg_x2{c}_run_mct) ',a,3e11.3,a,f9.6)"
    character(*),parameter :: subName = '(mrg_x2{c}_run_mct) '
!</list>

!<list> key="run_def" number=2
    integer  :: n, ki, ko, kir, kor
    integer  :: lsize
    real(r8) :: ifrac,ifracr
    real(r8) :: afrac,afracr
    real(r8) :: flux_epbalfact
    character(len=cl) :: flux_epbal
    type(seq_infodata_type),pointer :: infodata
    real(r8) :: frac_sum
    real(r8) :: avsdr, anidr, avsdf, anidf   ! albedos
    real(r8) :: fswabsv, fswabsi             ! sw

    !----- formats -----
    character(*),parameter :: F01 =   "('(mrg_x2{c}_run_mct) ',a,3e11.3,a,f9.6)"
    character(*),parameter :: subName = '(mrg_x2{c}_run_mct) '
!</list>

!<list> key="run_def" number=3
!wangty modify
#ifdef wrf
    integer  :: n, ki, kl, ko,k  ! indices, juanxiong he
#else
    integer  :: n, ki, kl, ko  ! indices
#endif
    real(r8) :: frac         ! temporary
    integer  :: lsize        ! temporary
!</list>

    !----------------------------------------------------------------------- 
    ! 
    ! Create input land state directly from atm output state
    !
!<list> key="prep" number=1
    call mct_avect_zero(x2{c}_{c})
    !
    ! Update surface fractions
    !
    ki=mct_aVect_indexRA(fractions_{c},"ifrac")
    kl=mct_aVect_indexRA(fractions_{c},"lfrac")
    ko=mct_aVect_indexRA(fractions_{c},"ofrac")
!</list>

#ifdef CPP_VECTOR
   usevector = .true.
#else
   usevector = .false.
#endif
!<list> key="cother" number=1
    call mct_aVect_copy(aVin=l2x_{c}, aVout=x2{c}_{c}, vector=usevector)
    call mct_aVect_copy(aVin=o2x_{c}, aVout=x2{c}_{c}, vector=usevector)
    call mct_aVect_copy(aVin=i2x_{c}, aVout=x2{c}_{c}, vector=usevector)
    call mct_aVect_copy(aVin=x{c}o_{c}, aVout=x2{c}_{c}, vector=usevector)
!</list>

!<list> key="cother" number=2
    call mct_aVect_copy(aVin=a2x_{c}, aVout=x2{c}_{c}, vector=usevector)
!</list>

!<list> key="cother" number=4
    call mct_aVect_copy(aVin=o2x_{c}, aVout=x2{c}_{c}, vector=usevector)
    call mct_aVect_copy(aVin=a2x_{c}, aVout=x2{c}_{c}, vector=usevector)
!</list>

!<list> key="cother" number=5
    call mct_aVect_copy(aVin=s2x_{c}, aVout=x2{c}_{c}, vector=usevector)
!</list>

!<list> key="cother" number=6
    call mct_aVect_copy(aVin=g2x_{c}, aVout=x2{c}_{c}, vector=usevector)
!</list>

!<list> key="run" number=1
    call seq_cdata_setptrs(cdata_{c1},infodata=infodata)
    ! Apply correction to precipitation of requested driver namelist

    call seq_infodata_GetData(infodata, flux_epbalfact = flux_epbalfact)

    ! Merge total snow and precip for ice input

    do i = 1,mct_aVect_lsize(x2{c1}_{c1})
       x2{c1}_{c1}%rAttr(index_x2{c1}_Faxa_rain,i) = a2x_{c1}%rAttr(index_a2x_Faxa_rainc,i) + &
                                                     a2x_{c1}%rAttr(index_a2x_Faxa_rainl,i)
       x2{c1}_{c1}%rAttr(index_x2{c1}_Faxa_snow,i) = a2x_{c1}%rAttr(index_a2x_Faxa_snowc,i) + &
                                                     a2x_{c1}%rAttr(index_a2x_Faxa_snowl,i)

       ! scale total precip and runoff by flux_epbalfact (TODO: note in cpl6
       ! this was always
       ! done over points where imask was > 0 - how does this translate here?)

       x2{c1}_{c1}%rAttr(index_x2{c1}_Faxa_rain,i) = x2{c1}_{c1}%rAttr(index_x2{c1}_Faxa_rain,i) * flux_epbalfact
       x2{c1}_{c1}%rAttr(index_x2{c1}_Faxa_snow,i) = x2{c1}_{c1}%rAttr(index_x2{c1}_Faxa_snow,i) * flux_epbalfact

    end do
!</list>

!<list> key="run" number=2
    call seq_cdata_setptrs(cdata_{c1}, infodata=infodata)
    call mct_aVect_zero(x2{c1}_{c1})

    call mct_aVect_copy(aVin=a2x_{c1}, aVout=x2{c1}_{c1}, vector=usevector)
    call mct_aVect_copy(aVin=i2x_{c1}, aVout=x2{c1}_{c1}, vector=usevector)
! tcx moved out to a separate accumulate
!!    call mct_aVect_copy(aVin=r2x_{c1}, aVout=x2{c1}_{c1}, vector=usevector)
    call mct_aVect_copy(aVin=xa{c1}_{c1}, aVout=x2{c1}_{c1}, vector=usevector)

    call seq_infodata_GetData(infodata, flux_epbalfact = flux_epbalfact)

    !
    ! Compute input ocn state (note that this only applies to non-land portion
    ! of gridcell)
    !
    ki  = mct_aVect_indexRa(fractions_{c1},"ifrac",perrWith=subName)
    ko  = mct_aVect_indexRa(fractions_{c1},"ofrac",perrWith=subName)
    kir = mct_aVect_indexRa(fractions_{c1},"ifrad",perrWith=subName)
    kor = mct_aVect_indexRa(fractions_{c1},"ofrad",perrWith=subName)
    lsize = mct_aVect_lsize(x2{c1}_{c1})
    do n = 1,lsize

       ifrac = fractions_{c1}%rAttr(ki,n)
       afrac = fractions_{c1}%rAttr(ko,n)
       frac_sum = ifrac + afrac
       if ((frac_sum) /= 0._r8) then
          ifrac = ifrac / (frac_sum)
          afrac = afrac / (frac_sum)
       endif

       ifracr = fractions_{c1}%rAttr(kir,n)
       afracr = fractions_{c1}%rAttr(kor,n)
       frac_sum = ifracr + afracr
       if ((frac_sum) /= 0._r8) then
          ifracr = ifracr / (frac_sum)
          afracr = afracr / (frac_sum)
       endif

       x2{c1}_{c1}%rAttr(index_x2{c1}_F{c1}xx_taux ,n) = xa{c1}_{c1}%rAttr(index_xa{c1}_Faox_taux ,n) * afrac + &
                                                      i2x_{c1}%rAttr(index_i2x_Fi{c1}i_taux ,n) * ifrac

       x2{c1}_{c1}%rAttr(index_x2{c1}_F{c1}xx_tauy ,n) = xa{c1}_{c1}%rAttr(index_xa{c1}_Faox_tauy ,n) * afrac + &
                                                      i2x_{c1}%rAttr(index_i2x_Fi{c1}i_tauy ,n) * ifrac

       ! --- was flux_solar:
       avsdr = xa{c1}_{c1}%rAttr(index_xa{c1}_So_avsdr,n)
       anidr = xa{c1}_{c1}%rAttr(index_xa{c1}_So_anidr,n)
       avsdf = xa{c1}_{c1}%rAttr(index_xa{c1}_So_avsdf,n)
       anidf = xa{c1}_{c1}%rAttr(index_xa{c1}_So_anidf,n)
       fswabsv  =  a2x_{c1}%rAttr(index_a2x_Faxa_swvdr,n) * (1.0_R8 - avsdr) &
                 + a2x_{c1}%rAttr(index_a2x_Faxa_swvdf,n) * (1.0_R8 - avsdf)
       fswabsi  =  a2x_{c1}%rAttr(index_a2x_Faxa_swndr,n) * (1.0_R8 - anidr) &
                 + a2x_{c1}%rAttr(index_a2x_Faxa_swndf,n) * (1.0_R8 - anidf)

       x2{c1}_{c1}%rAttr(index_x2{c1}_F{c1}xx_swnet,n) = (fswabsv + fswabsi) * afracr + &
                                                         i2x_{c1}%rAttr(index_i2x_Fi{c1}i_swpen,n) * ifrac

       x2{c1}_{c1}%rAttr(index_x2{c1}_F{c1}xx_lat  ,n) = xa{c1}_{c1}%rAttr(index_xa{c1}_Fa{c1}x_lat  ,n) * afrac

       x2{c1}_{c1}%rAttr(index_x2{c1}_F{c1}xx_sen  ,n) = xa{c1}_{c1}%rAttr(index_xa{c1}_Fa{c1}x_sen  ,n) * afrac

       x2{c1}_{c1}%rAttr(index_x2{c1}_F{c1}xx_evap ,n) = xa{c1}_{c1}%rAttr(index_xa{c1}_Fa{c1}x_evap ,n) * afrac

       x2{c1}_{c1}%rAttr(index_x2{c1}_F{c1}xx_lwup ,n) = xa{c1}_{c1}%rAttr(index_xa{c1}_Fa{c1}x_lwup ,n) * afrac

       x2{c1}_{c1}%rAttr(index_x2{c1}_F{c1}xx_lwdn ,n) = a2x_{c1}%rAttr(index_a2x_Faxa_lwdn ,n) * afrac

       x2{c1}_{c1}%rAttr(index_x2{c1}_F{c1}xx_snow ,n) = a2x_{c1}%rAttr(index_a2x_Faxa_snowc,n) * afrac + &
                                                         a2x_{c1}%rAttr(index_a2x_Faxa_snowl,n) * afrac

       x2{c1}_{c1}%rAttr(index_x2{c1}_F{c1}xx_rain ,n) = a2x_{c1}%rAttr(index_a2x_Faxa_rainc,n) * afrac + &
                                                         a2x_{c1}%rAttr(index_a2x_Faxa_rainl,n) * afrac

       x2{c1}_{c1}%rAttr(index_x2{c1}_F{c1}xx_melth,n) = i2x_{c1}%rAttr(index_i2x_Fi{c1}i_melth,n) * ifrac

       x2{c1}_{c1}%rAttr(index_x2{c1}_F{c1}xx_meltw,n) = i2x_{c1}%rAttr(index_i2x_Fi{c1}i_meltw,n) * ifrac

       x2{c1}_{c1}%rAttr(index_x2{c1}_F{c1}xx_salt ,n) = i2x_{c1}%rAttr(index_i2x_Fi{c1}i_salt ,n) * ifrac

       ! scale total precip and runoff by flux_epbalfact (TODO: note in cpl6 this was always
       ! done over points where imask was > 0 - how does this translate here?)

       x2{c1}_{c1}%rAttr(index_x2{c1}_F{c1}xx_rain ,n) = x2{c1}_{c1}%rAttr(index_x2{c1}_Foxx_rain ,n) * flux_epbalfact
       x2{c1}_{c1}%rAttr(index_x2{c1}_F{c1}xx_snow ,n) = x2{c1}_{c1}%rAttr(index_x2{c1}_Foxx_snow ,n) * flux_epbalfact
! this has been moved into r2xacc_rx
!       x2o_o%rAttr(index_x2o_Forr_roff ,n) = x2o_o%rAttr(index_x2o_Forr_roff
!       ,n) * flux_epbalfact
!       x2o_o%rAttr(index_x2o_Forr_ioff ,n) = x2o_o%rAttr(index_x2o_Forr_ioff
!       ,n) * flux_epbalfact

       x2{c1}_{c1}%rAttr(index_x2{c1}_F{c1}xx_prec ,n) = x2{c1}_{c1}%rAttr(index_x2{c1}_Foxx_rain ,n) + &
                                                      x2{c1}_{c1}%rAttr(index_x2{c1}_Foxx_snow ,n)

       x2{c1}_{c1}%rAttr(index_x2{c1}_F{c1}xx_bcphidry,n) = a2x_{c1}%rAttr(index_a2x_Faxa_bcphidry,n) * afrac
       x2{c1}_{c1}%rAttr(index_x2{c1}_F{c1}xx_bcphodry,n) = a2x_{c1}%rAttr(index_a2x_Faxa_bcphodry,n) * afrac
       x2{c1}_{c1}%rAttr(index_x2{c1}_F{c1}xx_bcphiwet,n) = a2x_{c1}%rAttr(index_a2x_Faxa_bcphiwet,n) * afrac
       x2{c1}_{c1}%rAttr(index_x2{c1}_F{c1}xx_ocphidry,n) = a2x_{c1}%rAttr(index_a2x_Faxa_ocphidry,n) * afrac
       x2{c1}_{c1}%rAttr(index_x2{c1}_F{c1}xx_ocphodry,n) = a2x_{c1}%rAttr(index_a2x_Faxa_ocphodry,n) * afrac
       x2{c1}_{c1}%rAttr(index_x2{c1}_F{c1}xx_ocphiwet,n) = a2x_{c1}%rAttr(index_a2x_Faxa_ocphiwet,n) * afrac
       x2{c1}_{c1}%rAttr(index_x2{c1}_F{c1}xx_dstwet1,n)  = a2x_{c1}%rAttr(index_a2x_Faxa_dstwet1,n)  * afrac
       x2{c1}_{c1}%rAttr(index_x2{c1}_F{c1}xx_dstwet2,n)  = a2x_{c1}%rAttr(index_a2x_Faxa_dstwet2,n)  * afrac
       x2{c1}_{c1}%rAttr(index_x2{c1}_F{c1}xx_dstwet3,n)  = a2x_{c1}%rAttr(index_a2x_Faxa_dstwet3,n)  * afrac
       x2{c1}_{c1}%rAttr(index_x2{c1}_F{c1}xx_dstwet4,n)  = a2x_{c1}%rAttr(index_a2x_Faxa_dstwet4,n)  * afrac
       x2{c1}_{c1}%rAttr(index_x2{c1}_F{c1}xx_dstdry1,n)  = a2x_{c1}%rAttr(index_a2x_Faxa_dstdry1,n)  * afrac
       x2{c1}_{c1}%rAttr(index_x2{c1}_F{c1}xx_dstdry2,n)  = a2x_{c1}%rAttr(index_a2x_Faxa_dstdry2,n)  * afrac
       x2{c1}_{c1}%rAttr(index_x2{c1}_F{c1}xx_dstdry3,n)  = a2x_{c1}%rAttr(index_a2x_Faxa_dstdry3,n)  * afrac
       x2{c1}_{c1}%rAttr(index_x2{c1}_F{c1}xx_dstdry4,n)  = a2x_{c1}%rAttr(index_a2x_Faxa_dstdry4,n)  * afrac

    end do

!</list>

!<list> key="run" number=3
    lsize = mct_avect_lsize(x2{c1}_{c1})
    do n = 1,lsize

       x2{c1}_{c1}%rAttr(index_x2{c1}_Sx_lfrac,n) = fractions_{c1}%Rattr(kl,n)
       x2{c1}_{c1}%rAttr(index_x2{c1}_Sx_ifrac,n) = fractions_{c1}%Rattr(ki,n)
       x2{c1}_{c1}%rAttr(index_x2{c1}_Sx_ofrac,n) = fractions_{c1}%Rattr(ko,n)

       frac = x2{c1}_{c1}%rAttr(index_x2{c1}_Sx_lfrac,n)
!wangty modify
#ifdef wrf
!--------------------------------------------------------------------------------------
! soil depth/height and soil temperature/moisture for wrf/cam coupling, added ! juanxiong he
!--------------------------------------------------------------------------------------
          do k=1,num_soil_layers
           x2{c1}_{c1}%rAttr(index_x2{c1}_Sx_soildepth(k),n)  = l2x_{c1}%rAttr(index_l2x_sl_soildepth(k),n)
           x2{c1}_{c1}%rAttr(index_x2{c1}_Sx_soilthick(k),n)  = l2x_{c1}%rAttr(index_l2x_sl_soilthick(k),n)
           x2{c1}_{c1}%rAttr(index_x2{c1}_Sx_soilt(k),n)  = l2x_{c1}%rAttr(index_l2x_sl_soilt(k),n)
           x2{c1}_{c1}%rAttr(index_x2{c1}_Sx_soilm(k),n)  = l2x_{c1}%rAttr(index_l2x_sl_soilm(k),n)
          enddo
!--------------------------------------------------------------------------------------
! soil depth/height and soil temperature/moisture for wrf/cam coupling, added ! juanxiong he
!--------------------------------------------------------------------------------------
#endif
       if (frac > 0._r8) then
          x2{c1}_{c1}%rAttr(index_x2{c1}_Sx_avsdr,n)  = x2{c1}_{c1}%rAttr(index_x2{c1}_Sx_avsdr,n) + &
                                               l2x_{c1}%rAttr(index_l2x_Sl_avsdr,n) * frac
          x2{c1}_{c1}%rAttr(index_x2{c1}_Sx_avsdf,n)  = x2{c1}_{c1}%rAttr(index_x2{c1}_Sx_avsdf,n) + &
                                               l2x_{c1}%rAttr(index_l2x_Sl_avsdf,n) * frac
          x2{c1}_{c1}%rAttr(index_x2{c1}_Sx_anidr,n)  = x2{c1}_{c1}%rAttr(index_x2{c1}_Sx_anidr,n) + &
                                               l2x_{c1}%rAttr(index_l2x_Sl_anidr,n) * frac
          x2{c1}_{c1}%rAttr(index_x2{c1}_Sx_anidf,n)  = x2{c1}_{c1}%rAttr(index_x2{c1}_Sx_anidf,n) + &
                                               l2x_{c1}%rAttr(index_l2x_Sl_anidf,n) * frac
          x2{c1}_{c1}%rAttr(index_x2{c1}_Sx_tref,n)   = x2{c1}_{c1}%rAttr(index_x2{c1}_Sx_tref,n) + &
                                               l2x_{c1}%rAttr(index_l2x_Sl_tref,n) * frac
          x2{c1}_{c1}%rAttr(index_x2{c1}_Sx_qref,n)   = x2{c1}_{c1}%rAttr(index_x2{c1}_Sx_qref,n) + &
                                               l2x_{c1}%rAttr(index_l2x_Sl_qref,n) * frac
          x2{c1}_{c1}%rAttr(index_x2{c1}_Sx_u10,n)    = x2{c1}_{c1}%rAttr(index_x2{c1}_Sx_u10,n) + &
                                               l2x_{c1}%rAttr(index_l2x_Sl_u10,n) * frac
          x2{c1}_{c1}%rAttr(index_x2{c1}_Faxx_evap,n) = x2{c1}_{c1}%rAttr(index_x2{c1}_Faxx_evap,n) + &
                                               l2x_{c1}%rAttr(index_l2x_Fall_evap,n) * frac
          x2{c1}_{c1}%rAttr(index_x2{c1}_Faxx_lwup,n) = x2{c1}_{c1}%rAttr(index_x2{c1}_Faxx_lwup,n) + &
                                               l2x_{c1}%rAttr(index_l2x_Fall_lwup,n) * frac
          x2{c1}_{c1}%rAttr(index_x2{c1}_Faxx_sen,n)  = x2{c1}_{c1}%rAttr(index_x2{c1}_Faxx_sen,n) + &
                                               l2x_{c1}%rAttr(index_l2x_Fall_sen,n) * frac
          x2{c1}_{c1}%rAttr(index_x2{c1}_Faxx_lat,n)  = x2{c1}_{c1}%rAttr(index_x2{c1}_Faxx_lat,n) + &
                                               l2x_{c1}%rAttr(index_l2x_Fall_lat,n) * frac
          x2{c1}_{c1}%rAttr(index_x2{c1}_Faxx_taux,n) = x2{c1}_{c1}%rAttr(index_x2{c1}_Faxx_taux,n) + &
                                               l2x_{c1}%rAttr(index_l2x_Fall_taux,n) * frac
          x2{c1}_{c1}%rAttr(index_x2{c1}_Faxx_tauy,n) = x2{c1}_{c1}%rAttr(index_x2{c1}_Faxx_tauy,n) + &
                                               l2x_{c1}%rAttr(index_l2x_Fall_tauy,n) * frac
          x2{c1}_{c1}%rAttr(index_x2{c1}_Sx_t,n)      = x2{c1}_{c1}%rAttr(index_x2{c1}_Sx_t,n) + &
                                               l2x_{c1}%rAttr(index_l2x_Sl_t,n) * frac
          !--- CO2 flux from lnd ---
          if ( (index_x2{c1}_Faxx_fco2_lnd /=0) .and. (index_l2x_Fall_fco2_lnd /=0)) then
             x2{c1}_{c1}%rAttr(index_x2{c1}_Faxx_fco2_lnd,n) = l2x_{c1}%rAttr(index_l2x_Fall_fco2_lnd,n) * frac
          end if

          if ( index_x2{c1}_Faxx_flxvoc1 /= 0 ) &
               x2{c1}_{c1}%rAttr(index_x2{c1}_Faxx_flxvoc1,n) = l2x_{c1}%rAttr(index_l2x_Fall_flxvoc1,n)* frac
          if ( index_x2{c1}_Faxx_flxvoc2 /= 0 ) &
               x2{c1}_{c1}%rAttr(index_x2{c1}_Faxx_flxvoc2,n) = l2x_{c1}%rAttr(index_l2x_Fall_flxvoc2,n)* frac
          if ( index_x2{c1}_Faxx_flxvoc3 /= 0 ) &
               x2{c1}_{c1}%rAttr(index_x2{c1}_Faxx_flxvoc3,n) = l2x_{c1}%rAttr(index_l2x_Fall_flxvoc3,n)* frac
          if ( index_x2{c1}_Faxx_flxvoc4 /= 0 ) &
               x2{c1}_{c1}%rAttr(index_x2{c1}_Faxx_flxvoc4,n) = l2x_{c1}%rAttr(index_l2x_Fall_flxvoc4,n)* frac
          if ( index_x2{c1}_Faxx_flxvoc5 /= 0 ) &
               x2{c1}_{c1}%rAttr(index_x2{c1}_Faxx_flxvoc5,n) = l2x_{c1}%rAttr(index_l2x_Fall_flxvoc5,n)* frac
       end if


       frac = x2{c1}_{c1}%rAttr(index_x2{c1}_Sx_ifrac,n)
       if (frac > 0._r8) then
          x2{c1}_{c1}%rAttr(index_x2{c1}_Sx_avsdr,n)  = x2{c1}_{c1}%rAttr(index_x2{c1}_Sx_avsdr,n) + &
                                                        i2x_{c1}%rAttr(index_i2x_Si_avsdr,n) * frac
          x2{c1}_{c1}%rAttr(index_x2{c1}_Sx_avsdf,n)  = x2{c1}_{c1}%rAttr(index_x2{c1}_Sx_avsdf,n) + &
                                                        i2x_{c1}%rAttr(index_i2x_Si_avsdf,n) * frac
          x2{c1}_{c1}%rAttr(index_x2{c1}_Sx_anidr,n)  = x2{c1}_{c1}%rAttr(index_x2{c1}_Sx_anidr,n) + &
                                                        i2x_{c1}%rAttr(index_i2x_Si_anidr,n) * frac
          x2{c1}_{c1}%rAttr(index_x2{c1}_Sx_anidf,n)  = x2{c1}_{c1}%rAttr(index_x2{c1}_Sx_anidf,n) + &
                                                        i2x_{c1}%rAttr(index_i2x_Si_anidf,n) * frac
          x2{c1}_{c1}%rAttr(index_x2{c1}_Sx_tref,n)   = x2{c1}_{c1}%rAttr(index_x2{c1}_Sx_tref,n) + &
                                                        i2x_{c1}%rAttr(index_i2x_Si_tref,n) * frac
          x2{c1}_{c1}%rAttr(index_x2{c1}_Sx_qref,n)   = x2{c1}_{c1}%rAttr(index_x2{c1}_Sx_qref,n) + &
                                                        i2x_{c1}%rAttr(index_i2x_Si_qref,n) * frac
          x2{c1}_{c1}%rAttr(index_x2{c1}_Sx_u10,n)    = x2{c1}_{c1}%rAttr(index_x2{c1}_Sx_u10,n) + &
                                                        i2x_{c1}%rAttr(index_i2x_Si_u10,n) * frac
          x2{c1}_{c1}%rAttr(index_x2{c1}_Faxx_evap,n) = x2{c1}_{c1}%rAttr(index_x2{c1}_Faxx_evap,n) + &
                                                        i2x_{c1}%rAttr(index_i2x_Faii_evap,n) * frac
          x2{c1}_{c1}%rAttr(index_x2{c1}_Faxx_lwup,n) = x2{c1}_{c1}%rAttr(index_x2{c1}_Faxx_lwup,n) + &
                                                        i2x_{c1}%rAttr(index_i2x_Faii_lwup,n) * frac
          x2{c1}_{c1}%rAttr(index_x2{c1}_Faxx_sen,n)  = x2{c1}_{c1}%rAttr(index_x2{c1}_Faxx_sen,n) + &
                                                        i2x_{c1}%rAttr(index_i2x_Faii_sen,n) * frac
          x2{c1}_{c1}%rAttr(index_x2{c1}_Faxx_lat,n)  = x2{c1}_{c1}%rAttr(index_x2{c1}_Faxx_lat,n) + &
                                                        i2x_{c1}%rAttr(index_i2x_Faii_lat,n) * frac
          x2{c1}_{c1}%rAttr(index_x2{c1}_Faxx_taux,n) = x2{c1}_{c1}%rAttr(index_x2{c1}_Faxx_taux,n) + &
                                                        i2x_{c1}%rAttr(index_i2x_Faii_taux,n) * frac
          x2{c1}_{c1}%rAttr(index_x2{c1}_Faxx_tauy,n) = x2{c1}_{c1}%rAttr(index_x2{c1}_Faxx_tauy,n) + &
                                                        i2x_{c1}%rAttr(index_i2x_Faii_tauy,n) * frac
          x2{c1}_{c1}%rAttr(index_x2{c1}_Sx_t,n)      = x2{c1}_{c1}%rAttr(index_x2{c1}_Sx_t,n) + &
                                                        i2x_{c1}%rAttr(index_i2x_Si_t,n) * frac
       end if

       frac = x2{c1}_{c1}%rAttr(index_x2{c1}_Sx_ofrac,n)
       if (frac > 0._r8) then
          x2{c1}_{c1}%rAttr(index_x2{c1}_Sx_avsdr,n)  = x2{c1}_{c1}%rAttr(index_x2{c1}_Sx_avsdr,n) + &
                                                        x{c1}o_{c1}%rAttr(index_x{c1}o_So_avsdr,n) * frac
          x2{c1}_{c1}%rAttr(index_x2{c1}_Sx_avsdf,n)  = x2{c1}_{c1}%rAttr(index_x2{c1}_Sx_avsdf,n) + &
                                                        x{c1}o_{c1}%rAttr(index_x{c1}o_So_avsdf,n) * frac
          x2{c1}_{c1}%rAttr(index_x2{c1}_Sx_anidr,n)  = x2{c1}_{c1}%rAttr(index_x2{c1}_Sx_anidr,n) + &
                                                        x{c1}o_{c1}%rAttr(index_x{c1}o_So_anidr,n) * frac
          x2{c1}_{c1}%rAttr(index_x2{c1}_Sx_anidf,n)  = x2{c1}_{c1}%rAttr(index_x2{c1}_Sx_anidf,n) + &
                                                        x{c1}o_{c1}%rAttr(index_x{c1}o_So_anidf,n) * frac
          x2{c1}_{c1}%rAttr(index_x2{c1}_Sx_tref,n)   = x2{c1}_{c1}%rAttr(index_x2{c1}_Sx_tref,n) + &
                                                        x{c1}o_{c1}%rAttr(index_x{c1}o_So_tref,n) * frac
          x2{c1}_{c1}%rAttr(index_x2{c1}_Sx_qref,n)   = x2{c1}_{c1}%rAttr(index_x2{c1}_Sx_qref,n) + &
                                                        x{c1}o_{c1}%rAttr(index_x{c1}o_So_qref,n) * frac
          x2{c1}_{c1}%rAttr(index_x2{c1}_Sx_u10,n)    = x2{c1}_{c1}%rAttr(index_x2{c1}_Sx_u10,n) + &
                                                        x{c1}o_{c1}%rAttr(index_x{c1}o_So_u10,n) * frac
          x2{c1}_{c1}%rAttr(index_x2{c1}_Faxx_evap,n) = x2{c1}_{c1}%rAttr(index_x2{c1}_Faxx_evap,n) + &
                                                        x{c1}o_{c1}%rAttr(index_x{c1}o_Faox_evap,n) * frac
          x2{c1}_{c1}%rAttr(index_x2{c1}_Faxx_lwup,n) = x2{c1}_{c1}%rAttr(index_x2{c1}_Faxx_lwup,n) + &
                                                        x{c1}o_{c1}%rAttr(index_x{c1}o_Faox_lwup,n) * frac
          x2{c1}_{c1}%rAttr(index_x2{c1}_Faxx_sen,n)  = x2{c1}_{c1}%rAttr(index_x2{c1}_Faxx_sen,n) + &
                                                        x{c1}o_{c1}%rAttr(index_x{c1}o_Faox_sen,n) * frac
          x2{c1}_{c1}%rAttr(index_x2{c1}_Faxx_lat,n)  = x2{c1}_{c1}%rAttr(index_x2{c1}_Faxx_lat,n) + &
                                                        x{c1}o_{c1}%rAttr(index_x{c1}o_Faox_lat,n) * frac
          x2{c1}_{c1}%rAttr(index_x2{c1}_Faxx_taux,n) = x2{c1}_{c1}%rAttr(index_x2{c1}_Faxx_taux,n) + &
                                                        x{c1}o_{c1}%rAttr(index_x{c1}o_Faox_taux,n) * frac
          x2{c1}_{c1}%rAttr(index_x2{c1}_Faxx_tauy,n) = x2{c1}_{c1}%rAttr(index_x2{c1}_Faxx_tauy,n) + &
                                                        x{c1}o_{c1}%rAttr(index_x{c1}o_Faox_tauy,n) * frac
          x2{c1}_{c1}%rAttr(index_x2{c1}_Sx_t,n)      = x2{c1}_{c1}%rAttr(index_x2{c1}_Sx_t,n) + &
                                                        o2x_a%rAttr(index_o2x_So_t,n) * frac
       end if
       !--- CO2 flux from ocn ---
       if ( (index_x2{c1}_Faxx_fco2_ocn /=0) .and. (index_o2x_Faoo_fco2 /=0)) then
          x2{c1}_{c1}%rAttr(index_x2{c1}_Faxx_fco2_ocn,n) = o2x_{c1}%rAttr(index_o2x_Faoo_fco2,n)
       end if

       !--- DMS flux from ocn ---
       if ( (index_x2{c1}_Faxx_fdms /=0) .and. (index_o2x_Faoo_fdms /=0)) then
          x2{c1}_{c1}%rAttr(index_x2{c1}_Faxx_fdms,n) = o2x_{c1}%rAttr(index_o2x_Faoo_fdms,n)
       end if

    end do

!</list>

  end subroutine mrg_x2{c}_run_mct
!
!===========================================================================================
!
  subroutine mrg_x2{c}_final_mct
    ! ******************
    ! Do nothing for now
    ! ******************
  end subroutine mrg_x2{c}_final_mct

end module mrg_x2{c}_mct

module seq_domain_mct

  use shr_kind_mod, only: R8=>shr_kind_r8, IN=>shr_kind_in
  use shr_kind_mod, only: CL=>shr_kind_cl
  use shr_sys_mod,  only: shr_sys_flush, shr_sys_abort
  use shr_mpi_mod,  only: shr_mpi_min, shr_mpi_max

  use mct_mod
  use seq_cdata_mod
  use seq_comm_mct
  use seq_infodata_mod

!{list} key="usemap"
!  use map_{ccc1}{ccc2}_mct

  implicit none
  save
  private ! except

!--------------------------------------------------------------------------
! Public interfaces
!--------------------------------------------------------------------------

  public :: seq_domain_check_mct
  public :: domain_areafactinit_mct

!--------------------------------------------------------------------------
! Public variables
!--------------------------------------------------------------------------

  real(R8), parameter :: eps_tiny   = 1.0e-16_R8 ! roundoff eps
  real(R8), parameter :: eps_big    = 1.0e+02_R8 ! big eps
  real(R8), parameter :: eps_frac_samegrid = 1.0e-14_R8 ! epsilon for fractions for samegrid

!--------------------------------------------------------------------------
! Private interfaces
!--------------------------------------------------------------------------

#ifdef CPP_VECTOR
  logical :: usevector = .true.
#else
  logical :: usevector = .false.
#endif
  
#ifdef SYSUNICOS
  logical :: usealltoall = .true.
#else
  logical :: usealltoall = .false.
#endif
  
  private :: seq_domain_check_grid_mct

#include <mpif.h>  

!================================================================================
contains
!================================================================================

!================================================================================
  subroutine seq_domain_check_mct(  &
!{list} key="listc"
!                                    cdata_{c}, &
                                 )
    !-----------------------------------------------------------
    !
    ! Arguments
    !
!{list} key="listc"
!    type(seq_cdata), intent(in) :: cdata_{c}
    !
    ! Local variables
    !
!{list} key="domainc"
!    type(mct_gGrid)  , pointer :: {ccc}dom_{c}   ! {ccc} domain
    !
!{list} key="listc"
!    type(mct_gsMap)  , pointer :: gsMap_{c}    ! {ccc} global seg map
    !
!{list} key="usemap"
!    type(mct_gGrid) :: {ccc2}dom_{c1}              ! {ccc2} domain info on {ccc1} decomp

    !
!{list} key="domainc"
!    integer(IN) :: mpicom_{c}                  ! {ccc} mpicom
    !
    type(seq_infodata_type), pointer :: infodata
    !
!{list} key="other_decomp"
!    real(R8), pointer :: frac{c}(:)            

!{list} key="other_decomp"
!    real(R8), pointer :: mask{c}(:)            

    !
    integer(IN) :: n, kl, ko, ki             ! indicies
    integer(IN) :: k1,k2,k3                  ! indicies
    !
!{list} key="present"
!    logical      :: {ccc}_present              ! {ccc} present flag

!{list} key="prognostic"
!    logical      :: {ccc1}{ccc2}_prognostic        ! {ccc1} {ccc2} prognostic flag

!{list} key="samegrid"
!    logical      :: samegrid_{c1}{c2}              ! {ccc1} {ccc2} grid same

    integer(IN)  :: rcode                    ! error status
!{list} key="domainc"
!    integer(IN)  :: {ccc}size                  ! local  size of {ccc}  grid

!{list} key="listc"
!    integer(IN)  :: g{ccc}size                 ! global size of {ccc}  grid

    integer(IN)  :: npts                     ! local size temporary
    integer(IN)  :: ier                      ! error code
    real(R8)     :: diff,dmaxo,dmaxi         ! difference tracker
    logical      :: iamroot                  ! local masterproc
    real(R8)     :: eps_frac                 ! epsilon for fractions
!{list} key="epsilon"
!    real(R8)     :: eps_{c1}{c2}mask               ! epsilon for masks
!{list} key="epsilon"
!    real(R8)     :: eps_{c1}{c2}grid               ! epsilon for grid coords
!{list} key="epsilon"
!    real(R8)     :: eps_{c1}{c2}area               ! epsilon for areas

    real(R8)     :: my_eps_frac              ! local eps_frac value
    real(r8)     :: rmin1,rmax1,rmin,rmax    ! local min max computation
    !
    real(R8),allocatable :: mask (:)         ! temporary real vector, domain mask
    !
    character(*),parameter :: F00 = "('(domain_check_mct) ',4a)"
    character(*),parameter :: F01 = "('(domain_check_mct) ',a,i6,a)"
    character(*),parameter :: F02 = "('(domain_check_mct) ',a,g23.15)"
    character(*),parameter :: F0R = "('(domain_check_mct) ',2A,2g23.15,A )"
    character(*),parameter :: subName = '(domain_check_mct) '
    !-----------------------------------------------------------

    call seq_comm_setptrs(CPLID,iamroot=iamroot)

    call seq_cdata_setptrs(cdata_a, infodata=infodata)
    call seq_infodata_GetData( infodata, &
!{list} key="samegrid"
!         samegrid_{c1}{c2}=samegrid_{c1}{c2}, &
         )
    call seq_infodata_GetData( infodata, &
!{list} key="present"
!         {ccc}_present={ccc}_present, &
!{list} key="prognostic"
!         {ccc1}{ccc2}_prognostic={ccc1}{ccc2}_prognostic, &
         )
    call seq_infodata_GetData( infodata, eps_frac=eps_frac, &
!{list} key="epsilon"
!        eps_{c1}mask=eps_{c1}{c2}mask, eps_{c1}grid=eps_{c1}{c2}grid, eps_{c1}area=eps_{c1}{c2}area, &
         )
    ! Get info

    call seq_cdata_setptrs(cdata_a, gsMap=gsMap_a, dom=atmdom_a, mpicom=mpicom_a)
    atmsize = mct_avect_lsize(atmdom_a%data)
    gatmsize = mct_gsMap_gsize(gsMap_a)
!<list> key="init" number=1
    if ({ccc2}_present) then
       call seq_cdata_setptrs(cdata_{c2}, gsMap=gsMap_{c2}, dom={ccc2}dom_{c2}, mpicom=mpicom_{c2})
       {ccc2}size = mct_avect_lsize({ccc2}dom_{c2}%data)
       g{ccc2}size = mct_gsMap_gsize(gsMap_{c2}) 
       if (samegrid_{c1}{c2} .and. g{ccc1}size /= g{ccc2}size) then
          write(logunit,*) subname,' error: global {ccc1}size = ',g{ccc1}size,' global {ccc2}size= ',g{ccc2}size
          call shr_sys_flush(logunit)
          call mct_die(subname,' {ccc1} and {ccc2} grid must have the same global size')
       end if
       call mct_gGrid_init(oGGrid={ccc2}dom_{c1}, iGGrid={ccc2}dom_{c2}, lsize={ccc1}size)
       call mct_aVect_zero({ccc2}dom_{c1}%data)
       call map_{ccc2}2{ccc1}_mct(cdata_{c2}, {ccc2}dom_{c2}%data, cdata_{c1}, {ccc2}dom_{c1}%data)
       allocate(mask{c2}({ccc1}size),stat=rcode)
       if(rcode /= 0) call mct_die(subName,'allocate mask{c2}')
       allocate(frac{c2}({ccc1}size),stat=rcode)
       if(rcode /= 0) call mct_die(subName,'allocate frac{c2}')
       call mct_aVect_exportRAttr({ccc2}dom_{c1}%data, 'mask', mask{c2}, {ccc1}size)
       call mct_aVect_exportRAttr({ccc2}dom_{c1}%data, 'frac', frac{c2}, {ccc1}size)
    endif
!</list>
!<list> key="init" number=1
    if ({ccc2}_present) then
       call seq_cdata_setptrs(cdata_{c2}, gsMap=gsMap_{c2}, dom={ccc2}dom_{c2}, mpicom=mpicom_{c2})
       {ccc2}size = mct_avect_lsize({ccc2}dom_{c2}%data)
       g{ccc2}size = mct_gsMap_gsize(gsMap_{c2})
       if (samegrid_{c1}{c2} .and. g{ccc1}size /= g{ccc2}size) then
          write(logunit,*) subname,' error: global {ccc1}size = ',g{ccc1}size,' global {ccc2}size= ',g{ccc2}size
          call shr_sys_flush(logunit)
          call mct_die(subname,' {ccc1} and {ccc2} grid must have the same global size')
       end if
       call mct_gGrid_init(oGGrid={ccc2}dom_{c1}, iGGrid={ccc2}dom_{c2}, lsize={ccc1}size)
       call mct_aVect_zero({ccc2}dom_{c1}%data)
       call map_{ccc2}2{ccc1}_mct(cdata_{c2}, {ccc2}dom_{c2}%data, cdata_{c1}, {ccc2}dom_{c1}%data)
       allocate(mask{c2}({ccc1}size),stat=rcode)
       if(rcode /= 0) call mct_die(subName,'allocate mask{c2}')
       allocate(frac{c2}({ccc1}size),stat=rcode)
       if(rcode /= 0) call mct_die(subName,'allocate frac{c2}')
       call mct_aVect_exportRAttr({ccc2}dom_{c1}%data, 'mask', mask{c2}, {ccc1}size)
       if (samegrid_{c1}{c2}) then
          call mct_aVect_exportRattr({ccc2}dom_{c1}%data, 'frac', frac{c2}, {ccc1}size)
       else
          call mct_aVect_exportRattr({ccc2}dom_{c1}%data, 'mask', frac{c2}, {ccc1}size)
       endif
    endif
!</list>

!<list> key="init" number=3
    if ({ccc1}_present .and. {ccc2}_present) then
       call seq_cdata_setptrs(cdata_{c2}, gsMap=gsMap_{c2}, dom={ccc2}dom_{c2}, mpicom=mpicom_{c2})
       {ccc2}size = mct_avect_lsize({ccc2}dom_{c2}%data)
       g{ccc2}size = mct_gsMap_gsize(gsMap_{c2})
       call seq_cdata_setptrs(cdata_{c1}, gsMap=gsMap_{c1}, dom={ccc1}dom_{c1}, mpicom=mpicom_{c1})
       {ccc1}size = mct_avect_lsize({ccc1}dom_{c1}%data)
       g{ccc1}size = mct_gsMap_gsize(gsMap_{c1})
       if (g{ccc2}size /= g{ccc2}size) then
          write(logunit,*) subname,' error: global {ccc2}size = ',g{ccc2}size,' global {ccc1}size= ',g{ccc1}size
          call shr_sys_flush(logunit)
          call mct_die(subname,' {ccc2} and {ccc1} grid must have the same global size')
       end if
       call mct_gGrid_init(oGGrid={ccc1}dom_{c2}, iGGrid={ccc1}dom_{c1}, lsize={ccc2}size)
       call mct_aVect_zero({ccc1}dom_{c2}%data)
       call map_{ccc1}2{ccc2}_mct(cdata_{c1}, {ccc1}dom_{c1}%data, cdata_{c2}, {ccc1}dom_{c2}%data)
       if (iamroot) write(logunit,F00) ' --- checking {ccc2}/{ccc1} domains ---'
       npts = {ccc2}size
       allocate(mask(npts),stat=rcode)
       if(rcode /= 0) call mct_die(subName,'allocate mask')
       call mct_aVect_getRAttr({ccc1}dom_{c2}%data,"mask",mask,rcode)
       where (mask < eps_axmask) mask = 0.0_r8
       call seq_domain_check_grid_mct({ccc2}dom_{c2}%data, {ccc1}dom_{c2}%data, 'mask', eps=eps_almask, mpicom=mpicom_{c2}, mask=mask)
       call seq_domain_check_grid_mct({ccc2}dom_{c2}%data, {ccc1}dom_{c2}%data, 'lat' , eps=eps_algrid, mpicom=mpicom_{c2}, mask=mask)
       call seq_domain_check_grid_mct({ccc2}dom_{c2}%data, {ccc1}dom_{c2}%data, 'lon' , eps=eps_algrid, mpicom=mpicom_{c2}, mask=mask)
       call seq_domain_check_grid_mct({ccc2}dom_{c2}%data, {ccc1}dom_{c2}%data, 'area', eps=eps_alarea, mpicom=mpicom_{c2}, mask=mask)
       deallocate(mask,stat=rcode)
       if(rcode /= 0) call mct_die(subName,'deallocate mask')
    endif
!</list>

!<list> key="init" number=4
    if ({ccc1}_present .and. {ccc2}_present) then
       if (g{ccc2}size /= g{ccc1}size) then
          write(logunit,*) subname,' error: global {ccc2}size = ',g{ccc2}size,' global {ccc1}size= ',g{ccc1}size
          call shr_sys_flush(logunit)
          call mct_die(subname,' {ccc2} and {ccc1} grid must have the same global
size')
       endif
       call mct_gGrid_init(oGGrid={ccc1}dom_{c2}, iGGrid={ccc1}dom_{c1}, lsize={ccc2}size)
       call mct_aVect_zero({ccc1}dom_{c2}%data)
       call map_{ccc1}2{ccc2}_mct(cdata_{c1}, {ccc1}dom_{c1}%data, cdata_{c2}, {ccc1}dom_{c2}%data)
    end if
!</list>

!<list> key="init" number=5
    if ({ccc1}_present .and. {ccc2}{ccc1}_prognostic .and. samegrid_{c1}{c2}) then
       call seq_cdata_setptrs(cdata_{c1}, gsMap=gsMap_{c1})
       g{ccc1}size = mct_gsMap_gsize(gsMap_{c1})
       if (g{ccc2}size /= g{ccc1}size) then
          write(logunit,*) subname,' error: global {ccc2}size = ',g{ccc2}size,' global {ccc1}size= ',g{ccc1}size
          call shr_sys_flush(logunit)
          call mct_die(subname,' {ccc2} and {ccc1} grid must have the same global size')
       endif
    end if
!</list>
!<list> key="init" number=4
    !------------------------------------------------------------------------------
    ! Check {ccc1}/{ccc2} grid consistency
    !------------------------------------------------------------------------------
     if ({ccc2}_present .and. {ccc1}_present) then
!    if (samegrid_oi) then       ! doesn't yet exist

       npts = {ccc2}size
       allocate(mask(npts),stat=rcode)
       if(rcode /= 0) call mct_die(subName,'allocate mask')

       if (iamroot) write(logunit,F00) ' --- checking {ccc2}/{ccc1} domains ---'
       call seq_domain_check_grid_mct({ccc2}dom_{c2}%data, {ccc1}dom_{c2}%data,'mask', eps=eps_{c2}{c1}grid, mpicom=mpicom_{c2})
       call mct_aVect_getRAttr({ccc2}dom_{c2}%data,"mask",mask,rcode)
       where (mask < eps_{c2}{c1}mask) mask = 0.0_r8

       call seq_domain_check_grid_mct({ccc2}dom_{c2}%data, {ccc1}dom_{c2}%data,'lat' , eps=eps_{c2}{c1}grid, mpicom=mpicom_{c2}, mask=mask)
       call seq_domain_check_grid_mct({ccc2}dom_{c2}%data, {ccc1}dom_{c2}%data,'lon' , eps=eps_{c2}{c1}grid, mpicom=mpicom_{c2}, mask=mask)
       call seq_domain_check_grid_mct({ccc2}dom_{c2}%data, {ccc1}dom_{c2}%data,'area', eps=eps_{c2}{c1}grid, mpicom=mpicom_{c2}, mask=mask)

       deallocate(mask,stat=rcode)
       if(rcode /= 0) call mct_die(subName,'deallocate mask')

!    endif
     endif
!</list>

!<list> key="init" number=1
    !------------------------------------------------------------------------------
    ! Check {ccc1}/{ccc2} grid consistency
    !------------------------------------------------------------------------------

    if ({ccc2}_present .and. samegrid_{c1}{c2}) then
       if (iamroot) write(logunit,F00) ' --- checking {ccc1}/{ccc2} domains ---'
       call seq_domain_check_grid_mct({ccc1}dom_{c1}%data, {ccc2}dom_{c1}%data, 'lat' , eps=eps_{c1}{c2}grid, mpicom=mpicom_{c1}, mask=mask{c2})
       call seq_domain_check_grid_mct({ccc1}dom_{c1}%data, {ccc2}dom_{c1}%data, 'lon' , eps=eps_{c1}{c2}grid, mpicom=mpicom_{c1}, mask=mask{c2})
       call seq_domain_check_grid_mct({ccc1}dom_{c1}%data, {ccc2}dom_{c1}%data, 'area', eps=eps_{c1}{c2}grid, mpicom=mpicom_{c1}, mask=mask{c2})

    endif
!</list>

!<list> key="init" number=2
    !------------------------------------------------------------------------------
    ! Check {ccc1}/{ccc2} grid consistency (if samegrid)
    !------------------------------------------------------------------------------

    if ({ccc2}_present .and. samegrid_{c1}o) then
       if (iamroot) write(logunit,F00) ' --- checking {ccc1}/{ccc2} domains ---'
       call seq_domain_check_grid_mct({ccc1}dom_{c1}%data, {ccc2}dom_{c1}%data, 'lat' , eps=eps_algrid, mpicom=mpicom_{c1}, mask=mask{c2})
       call seq_domain_check_grid_mct({ccc1}dom_{c1}%data, {ccc2}dom_{c1}%data, 'lon' , eps=eps_algrid, mpicom=mpicom_{c1}, mask=mask{c2})
       call seq_domain_check_grid_mct({ccc1}dom_{c1}%data, {ccc2}dom_{c1}%data, 'area', eps=eps_algrid, mpicom=mpicom_{c1}, mask=mask{c2})
    endif
!</list>

    !------------------------------------------------------------------------------
    ! Check consistency of land fraction with ocean mask on grid
    !------------------------------------------------------------------------------

    my_eps_frac = eps_frac
    if (samegrid_ao) my_eps_frac = eps_frac_samegrid
    if (.not. samegrid_al) my_eps_frac = eps_big

    if (iamroot) write(logunit,F00) ' --- checking fractions in domains ---'
    dmaxi = 0.0_R8
    dmaxo = 0.0_R8
    do n = 1,atmsize
       if (lnd_present .and. ice_present) then
          diff = abs(1._r8 - fracl(n) - fraci(n))
          dmaxi = max(diff,dmaxi)
          if (diff > my_eps_frac) then
             write(logunit,*)'inconsistency between land fraction and sea ice fraction'
             write(logunit,*)'n= ',n,' fracl= ',fracl(n),' fraci= ',fraci(n),' sum= ',fracl(n)+fraci(n)
             call shr_sys_flush(logunit)
             call mct_die(subname,'inconsistency between land fraction and sea ice fraction')
          end if
          if ((1._r8-fraci(n)) > eps_frac .and. fracl(n) < eps_tiny) then
             write(logunit,*)'inconsistency between land mask and sea ice mask'
             write(logunit,*)'n= ',n,' fracl= ',fracl(n),' fraci= ',fraci(n)
             call shr_sys_flush(logunit)
             call mct_die(subname,' inconsistency between land mask and sea ice mask')
          end if
       endif
       if (lnd_present .and. ocn_present) then
          diff = abs(1._r8 - fracl(n) - fraco(n))
          dmaxo = max(diff,dmaxo)
          if (diff > my_eps_frac) then
             write(logunit,*)'inconsistency between land fraction and ocn land fraction'
             write(logunit,*)'n= ',n,' fracl= ',fracl(n),' fraco= ',fraco(n),' sum= ',fracl(n)+fraco(n)
             call shr_sys_flush(logunit)
             call mct_die(subname,' inconsistency between land fraction and ocn land fraction')
          end if
          if ((1._r8-fraco(n)) > eps_frac .and. fracl(n) < eps_tiny) then
             write(logunit,*)'inconsistency between land mask and ocn land mask'
             write(logunit,*)'n= ',n,' fracl= ',fracl(n),' fraco= ',fraco(n)
             call shr_sys_flush(logunit)
             call mct_die(subname,' inconsistency between land mask and ocn land mask')
          end if
       endif
    end do 
    if (iamroot) then
       write(logunit,F02) ' maximum           difference for ofrac sum ',dmaxo
       write(logunit,F02) ' maximum           difference for ifrac sum ',dmaxi
       write(logunit,F02) ' maximum allowable difference for  frac sum ',my_eps_frac
       write(logunit,F02) ' maximum allowable tolerance for valid frac ',eps_frac
       call shr_sys_flush(logunit)
    endif

    !------------------------------------------------------------------------------
    ! Set atm and lnd ascale
    !------------------------------------------------------------------------------
    if (lnd_present .and. ocn_present) then
       k1 = mct_aVect_indexRa(atmdom_a%data,"ascale",perrWith='domain_check ascale')
       do k2 = 1,atmsize
          if (fracl(k2) /= 0.0_r8) then
             atmdom_a%data%rAttr(k1,k2) = (1.0_r8-fraco(k2))/(fracl(k2))
          else
             atmdom_a%data%rAttr(k1,k2) = 0.0
          endif
       enddo
       call map_atm2lnd_mct(cdata_a, atmdom_a%data, cdata_l, lnddom_l%data, fluxlist='ascale')
    elseif (lnd_present .and. ice_present) then
       k1 = mct_aVect_indexRa(atmdom_a%data,"ascale",perrWith='domain_check ascale')
       do k2 = 1,atmsize
          if (fracl(k2) /= 0.0_r8) then
             atmdom_a%data%rAttr(k1,k2) = (1.0_r8-fraci(k2))/(fracl(k2))
          else
             atmdom_a%data%rAttr(k1,k2) = 0.0
          endif
       enddo
       call map_atm2lnd_mct(cdata_a, atmdom_a%data, cdata_l, lnddom_l%data, fluxlist='ascale')
    endif

    if (lnd_present .and. (ocn_present .or. ice_present)) then
       k1 = mct_aVect_indexRa(atmdom_a%data,"ascale",perrWith='domain_check atm ascale')
       rmin1 = minval(atmdom_a%data%rAttr(k1,:))
       rmax1 = maxval(atmdom_a%data%rAttr(k1,:))
       call shr_mpi_min(rmin1,rmin,mpicom_a)
       call shr_mpi_max(rmax1,rmax,mpicom_a)
       if (iamroot) write(logunit,F0R) trim(subname),' : min/max ascale ',rmin,rmax,' atmdom_a'

       k1 = mct_aVect_indexRa(lnddom_l%data,"ascale",perrWith='domain_check lnd ascale')
       rmin1 = minval(lnddom_l%data%rAttr(k1,:))
       rmax1 = maxval(lnddom_l%data%rAttr(k1,:))
       call shr_mpi_min(rmin1,rmin,mpicom_l)
       call shr_mpi_max(rmax1,rmax,mpicom_l)
       if (iamroot) write(logunit,F0R) trim(subname),' : min/max ascale ',rmin,rmax,' lnddom_l'
    endif

    !------------------------------------------------------------------------------
    ! Clean up allocated memory
    !------------------------------------------------------------------------------

    if (lnd_present) then
       deallocate(fracl,stat=rcode)
       if(rcode /= 0) call mct_die(subName,'deallocate fracl')
       deallocate(maskl,stat=rcode)
       if(rcode /= 0) call mct_die(subName,'deallocate maskl')
       call mct_gGrid_clean(lnddom_a, rcode)
       if(rcode /= 0) call mct_die(subName,'clean lnddom_a')
    endif

    if (ocn_present) then
       deallocate(fraco,stat=rcode)
       if(rcode /= 0) call mct_die(subName,'deallocate fraco')
       deallocate(masko,stat=rcode)
       if(rcode /= 0) call mct_die(subName,'deallocate masko')
       call mct_gGrid_clean(ocndom_a, rcode)
       if(rcode /= 0) call mct_die(subName,'clean ocndom_a')
    endif

    if (ice_present) then
       deallocate(fraci,stat=rcode)
       if(rcode /= 0) call mct_die(subName,'deallocate fraci')
       deallocate(maski,stat=rcode)
       if(rcode /= 0) call mct_die(subName,'deallocate maski')
       call mct_gGrid_clean(icedom_a, rcode)
       if(rcode /= 0) call mct_die(subName,'clean icedom_o')
    endif

    if (ocn_present .and. ice_present) then
       call mct_gGrid_clean(icedom_o, rcode)
       if(rcode /= 0) call mct_die(subName,'clean icedom_o')
    endif

  end subroutine seq_domain_check_mct

!===============================================================================
  
  subroutine seq_domain_check_grid_mct(dom1, dom2, attr, eps, mpicom, mask)
   
    !-----------------------------------------------------------

    ! Arguments

    type(mct_aVect) , intent(in) :: dom1
    type(mct_aVect) , intent(in) :: dom2
    character(len=*), intent(in) :: attr   ! grid attribute to compare
    real(R8)        , intent(in) :: eps    ! error condition for compare
    integer(IN)     , intent(in) :: mpicom
    real(R8)        , intent(in), optional :: mask(:)

    ! Local variables

    integer(in)       :: n,ndiff            ! indices
    integer(in)       :: npts1,npts2,npts   ! counters
    integer(in)       :: rcode              ! error code
    real(R8)          :: diff,max_diff      ! temporaries
    real(R8)          :: tot_diff           ! maximum diff across all pes
    integer(IN)       :: ier                ! error code
    real(R8), pointer :: data1(:)           ! temporaries
    real(R8), pointer :: data2(:)           ! temporaries
    real(R8), pointer :: lmask(:)           ! temporaries
    logical           :: iamroot            ! local masterproc

    character(*),parameter :: F00 = "('(domain_check_grid_mct) ',4a)"
    character(*),parameter :: F01 = "('(domain_check_grid_mct) ',a,i12,a)"
    character(*),parameter :: F02 = "('(domain_check_grid_mct) ',2a,g23.15)"
    character(*),parameter :: F0R = "('(domain_check_grid_mct) ',2A,2g23.15,A )"
    character(*),parameter :: subName = '(domain_check_grid_mct) '
    !-----------------------------------------------------------

    call seq_comm_setptrs(CPLID,iamroot=iamroot)

   npts1 = mct_aVect_lsize(dom1)
   npts2 = mct_aVect_lsize(dom2)
   npts  = npts1
   
   if (npts1 == npts2) then
      if (iamroot) write(logunit,F01) " the domain size is = ", npts
   else
      write(logunit,*) trim(subname)," domain size #1 = ", npts1
      write(logunit,*) trim(subname)," domain size #2 = ", npts2
      write(logunit,*) trim(subname)," ERROR: domain size mis-match"
      call mct_die(subName,"ERROR: domain size mis-match")
   end if
   
   allocate(data1(npts),stat=rcode)
   if(rcode /= 0) call mct_die(subName,'allocate data1')
   allocate(data2(npts),stat=rcode)
   if(rcode /= 0) call mct_die(subName,'allocate data2')
   allocate(lmask(npts),stat=rcode)
   if(rcode /= 0) call mct_die(subName,'allocate lmask')

   call mct_aVect_exportRAttr(dom1, trim(attr), data1, npts)
   call mct_aVect_exportRAttr(dom2, trim(attr), data2, npts)
   lmask = 1.0_r8
   if (present(mask)) then
      if (size(mask) /= npts) then
         call mct_die(subName,"ERROR: mask size mis-match")
      endif
      lmask = mask
   endif

   ! --- adjust lons to address wraparound issues, we're assuming degree here! ---

   if (trim(attr) == "lon") then
      do n = 1,npts
         if (data2(n) > data1(n)) then
            do while ( (data1(n)+360.0_R8) < (data2(n)+180.0_R8) ) ! longitude is periodic
               data1(n) = data1(n) + 360.0_R8
            end do
         else
            do while ( (data2(n)+360.0_R8) < (data1(n)+180.0_R8) ) ! longitude is periodic
               data2(n) = data2(n) + 360.0_R8
            end do
         endif
      enddo
   endif

   ! Only check consistency where mask is greater than zero, if mask is present

   max_diff = 0.0_R8
   ndiff = 0
   do n=1,npts
      if (lmask(n) > eps_tiny) then
         diff = abs(data1(n)-data2(n))
         max_diff = max(max_diff,diff)
         if (diff > eps) then
!debug            write(logunit,*)'n= ',n,' data1= ',data1(n),' data2= ',data2(n),' diff= ',diff, ' eps= ',eps
            ndiff = ndiff + 1
         endif
      end if
   end do

   call mpi_reduce(max_diff,tot_diff,1,MPI_REAL8,MPI_MAX,0,mpicom,ier)
   if (iamroot) then
      write(logunit,F02) " maximum           difference for ",trim(attr),tot_diff
      write(logunit,F02) " maximum allowable difference for ",trim(attr),eps
      call shr_sys_flush(logunit)
   endif
   call mpi_barrier(mpicom,ier)

   if (ndiff > 0) then
      write(logunit,*) trim(subname)," ERROR: incompatible domain grid coordinates"
      call shr_sys_flush(logunit)
      call mct_die(subName,"incompatible domain grid coordinates")
   endif
   
   deallocate(data1,stat=rcode)
   if(rcode /= 0) call mct_die(subName,'deallocate data1')
   deallocate(data2,stat=rcode)
   if(rcode /= 0) call mct_die(subName,'deallocate data2')
   deallocate(lmask,stat=rcode)
   if(rcode /= 0) call mct_die(subName,'deallocate lmask')

 end subroutine seq_domain_check_grid_mct

!===============================================================================

 subroutine domain_areafactinit_mct( cdata, mdl2drv, drv2mdl, comment)
    !-----------------------------------------------------------
    !
    ! Arguments
    !
    type(seq_cdata) , intent(inout) :: cdata
    real(R8),pointer                :: mdl2drv(:)
    real(R8),pointer                :: drv2mdl(:)
    character(len=*),optional,intent(in) :: comment
    !
    ! Local variables
    !
    type(seq_infodata_type),pointer :: infodata
    type(mct_gGrid),pointer:: domain
    integer                :: ID
    integer                :: mpicom
    logical                :: iamroot
    logical                :: samegrid_ao, samegrid_al
    integer                :: j1,j2,m1,n,rcode
    integer                :: gridsize,m2dsize,d2msize
    real(r8)               :: rmin1,rmax1,rmin,rmax
    real(r8)               :: rmask,rarea,raream
    character(cl)          :: lcomment
    character(len=*),parameter :: subName = '(domain_areafactinit_mct) '
    character(len=*),parameter :: F0R = "(2A,2g23.15,A )"
    !
    !-----------------------------------------------------------

    lcomment = ''
    if (present(comment)) lcomment = comment

    call seq_cdata_setptrs(cdata, ID=ID, dom=domain, infodata=infodata)
    call seq_comm_setptrs(ID, mpicom=mpicom, iamroot=iamroot) 
    call seq_infodata_GetData( infodata, samegrid_ao=samegrid_ao, samegrid_al=samegrid_al)

    ! get sizes

    gridsize = mct_gGrid_lsize(domain)
    allocate(drv2mdl(gridsize),mdl2drv(gridsize),stat=rcode)
    if(rcode /= 0) call mct_die(subname,'allocate area correction factors')

    j1 = mct_gGrid_indexRA(domain,"area"    ,dieWith=subName)
    j2 = mct_gGrid_indexRA(domain,"aream"   ,dieWith=subName)
    m1 = mct_gGrid_indexRA(domain,"mask"    ,dieWith=subName)

    mdl2drv(:)=1.0_R8
    drv2mdl(:)=1.0_R8

    if (samegrid_ao .and. samegrid_al) then
        ! default 1.0
    else
       do n=1,gridsize
          rmask  = domain%data%rAttr(m1,n)
          rarea  = domain%data%rAttr(j1,n)
          raream = domain%data%rAttr(j2,n)
          if ( abs(rmask) >= 1.0e-06) then
             if (rarea * raream /= 0.0_r8) then
                mdl2drv(n) = rarea/raream
                drv2mdl(n) = 1.0_R8/mdl2drv(n)
                !if (mdl2drv(n) > 10.0 .or. mdl2drv(n) < 0.1) then
                !   write(logunit,*) trim(subname),' WARNING area,aream= ', &
                !      domain%data%rAttr(j1,n),domain%data%rAttr(j2,n),' in ',n,gridsize
                !endif
             else
                write(logunit,*) trim(subname),' ERROR area,aream= ', &
                   rarea,raream,' in ',n,gridsize
                call shr_sys_flush(logunit)
                call shr_sys_abort()
             endif
          endif
       enddo
    end if
       
    rmin1 = minval(mdl2drv)
    rmax1 = maxval(mdl2drv)
    call shr_mpi_min(rmin1,rmin,mpicom)
    call shr_mpi_max(rmax1,rmax,mpicom)
    if (iamroot) write(logunit,F0R) trim(subname),' : min/max mdl2drv ',rmin,rmax,trim(lcomment)

    rmin1 = minval(drv2mdl)
    rmax1 = maxval(drv2mdl)
    call shr_mpi_min(rmin1,rmin,mpicom)
    call shr_mpi_max(rmax1,rmax,mpicom)
    if (iamroot) write(logunit,F0R) trim(subname),' : min/max drv2mdl ',rmin,rmax,trim(lcomment)

 end subroutine domain_areafactinit_mct

!===============================================================================

end module seq_domain_mct




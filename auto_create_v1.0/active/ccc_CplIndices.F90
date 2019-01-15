module {ccc}_CplIndices
  
  use seq_flds_mod
  use mct_mod
!xmlinsert(list[module_use])

  implicit none

  SAVE
  public                               ! By default make data private
!xmlinsert(list[var_model_to_cpl])
!xmlinsert(list[var_cpl_to_model])

!{list} key="var_model_to_cpl"

!{list} key="var_cpl_to_model"

contains

!{list} key="indices_set" number=1
!  subroutine {ccc}_CplIndicesSet( )

!{list} key="indices_set" number=2
!  subroutine {ccc}_cpl_indices_set( )

!xmlinsert(list[indices_set_use])
    type(mct_aVect) :: {c}2x      ! temporary
    type(mct_aVect) :: x2{c}      ! temporary
!xmlinsert(list[indices_set_var])

    ! Determine attribute vector indices

    ! create temporary attribute vectors
    call mct_aVect_init(x2{c}, rList=seq_flds_x2{c}_fields, lsize=1)
    call mct_aVect_init({c}2x, rList=seq_flds_{c}2x_fields, lsize=1)

!xmlinsert(list[other_init])

!{list} key="init_var_x2c"
!xmlinsert(list[init_var_x2c])

!{list} key="init_var_c2x"
!xmlinsert(list[init_var_c2x])

    call mct_aVect_clean(x2{c})
    call mct_aVect_clean({c}2x)
    !call mct_aVect_clean(r2x_o) !LPF 20121219

!xmlinsert(list[other_clean])

  end subroutine {ccc}_CplIndicesSet

end module {ccc}_CplIndices


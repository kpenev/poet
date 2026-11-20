! ***********************************************************************
!
!   Merged run_star_extras: MIST2 base (r7503-native) with POET history
!   columns (rcore/mcore/core_inertia/rrad/mrad/rad_inertia/env_inertia)
!   and POET wind subroutines (Reimers_then_Blocker / Reimers_then_VW).
!
! ***********************************************************************
module run_star_extras

  use star_lib
  use star_def
  use const_def
  use crlibm_lib

  implicit none

  real(dp) :: original_diffusion_dt_limit
  logical :: rot_set_check = .false.
  logical :: TP_AGB_check = .false.
  logical :: late_AGB_check = .false.
  logical :: post_AGB_check = .false.
  logical :: pre_WD_check = .false.

contains

  subroutine extras_controls(id, ierr)
    integer, intent(in) :: id
    integer, intent(out) :: ierr
    type (star_info), pointer :: s
    ierr = 0
    call star_ptr(id, s, ierr)
    if (ierr /= 0) return

    s% other_am_mixing => TSF
    s% other_wind => Reimers_then_Blocker

    s% extras_startup => extras_startup
    s% extras_check_model => extras_check_model
    s% extras_finish_step => extras_finish_step
    s% extras_after_evolve => extras_after_evolve
    s% how_many_extra_history_columns => how_many_extra_history_columns
    s% data_for_extra_history_columns => data_for_extra_history_columns
    s% how_many_extra_profile_columns => how_many_extra_profile_columns
    s% data_for_extra_profile_columns => data_for_extra_profile_columns

    s% job% warn_run_star_extras = .false.

    original_diffusion_dt_limit = s% diffusion_dt_limit

  end subroutine extras_controls


  integer function extras_startup(id, restart, ierr)
    integer, intent(in) :: id
    logical, intent(in) :: restart
    integer, intent(out) :: ierr
    type (star_info), pointer :: s
    integer :: j, cid
    real(dp) :: frac, vct30, vct100, scale_factor
    ierr = 0
    call star_ptr(id, s, ierr)
    if (ierr /= 0) return
    extras_startup = 0
    if (.not. restart) then
       call alloc_extra_info(s)
    else
       call unpack_extra_info(s)
    end if

    if (.not. s% job% extras_lpar(1)) rot_set_check = .true.

    if(s% initial_mass < 10.0d0 .and. s% initial_mass >= 0.6d0)then
       TP_AGB_check=.true.
    endif

    if(s% initial_mass < 7.0d0) then
       if(s% initial_mass > 4.0d0) then
          scale_factor = pow2((7d0 - s% initial_mass)/(7d0 - 4d0))
       else
          scale_factor = 1d0
       end if

       s% am_nu_ST_factor = 0.0d0
       s% use_other_am_mixing = .true.
       s% am_time_average = .true.
       s% premix_omega = .true.
       s% recalc_mixing_info_each_substep = .true.
       s% am_nu_factor = 1.0d0
       s% am_nu_non_rotation_factor = 1.0d0
       s% angsml = 0.0d0
       s% am_nu_visc_factor = 0.333d0 + (1d0 - scale_factor)*(1d0 - 0.333d0)

       print *, "Applying Fuller et al Tayler-Spruit with scale factor", scale_factor
       print *, "am_nu_visc_factor", s% am_nu_visc_factor
    endif

    s% overshoot_f_above_burn_h_core  = 0.016_dp
    s% overshoot_f0_above_burn_h_core = 8.0d-3

    ! extras_rpar(4) = [Fe/H]; set correctly by inlist patching
    s% overshoot_f_below_nonburn_shell  = f_ov_below_nonburn(s% job% extras_rpar(4), s% initial_mass)
    s% overshoot_f0_below_nonburn_shell = 0.5d0 * s% overshoot_f_below_nonburn_shell

  end function extras_startup


  function f_ov_below_nonburn(feh, minit) result(f_ov)
    real(dp), intent(in) :: feh, minit
    real(dp) :: f_ov
    real(dp), parameter :: max_f = 8.0d-2
    real(dp), parameter :: min_f = 1.0d-2
    if (minit < 2.0_dp) then
      f_ov = 0.016_dp - 0.027_dp*feh
      f_ov = min(max(f_ov, min_f), max_f)
    else
      f_ov = 0.016_dp
    endif
  end function f_ov_below_nonburn


  function f_ov_fcn_of_mass(m) result(f_ov)
    real(dp), intent(in) :: m
    real(dp) :: f_ov, frac
    real(dp), parameter :: f1 = 0.016_dp, f2 = 0.0415_dp
    if(m < 4.0d0) then
       frac = 0.0d0
    elseif(m > 8.0d0) then
       frac = 1.0d0
    else
       frac = 0.5_dp * (1.0_dp - cos( pi * 0.25_dp*(m - 4._dp) ))
    endif
    f_ov = f1 + (f2-f1)*frac
  end function f_ov_fcn_of_mass


  integer function extras_check_model(id, id_extra)
    integer, intent(in) :: id, id_extra
    integer :: ierr
    type (star_info), pointer :: s
    ierr = 0
    call star_ptr(id, s, ierr)
    if (ierr /= 0) return
    extras_check_model = keep_going
  end function extras_check_model


  ! -----------------------------------------------------------------------
  ! Extra history columns — POET layout (9 columns):
  !   1 rcore        radius at convective-core / radiative-shell boundary
  !   2 mcore        mass   at convective-core / radiative-shell boundary
  !   3 core_inertia moment of inertia of convective core
  !   4 rrad         radius at radiative-shell / convective-envelope boundary
  !   5 mrad         mass   at radiative-shell / convective-envelope boundary
  !   6 rad_inertia  moment of inertia of radiative shell
  !   7 env_inertia  moment of inertia of convective envelope
  !   8 core_outer_i diagnostic: zone index of core outer boundary
  !   9 env_bottom_i diagnostic: zone index of envelope bottom
  ! -----------------------------------------------------------------------

  integer function how_many_extra_history_columns(id, id_extra)
    integer, intent(in) :: id, id_extra
    integer :: ierr
    type (star_info), pointer :: s
    ierr = 0
    call star_ptr(id, s, ierr)
    if (ierr /= 0) return
    how_many_extra_history_columns = 9
  end function how_many_extra_history_columns


  ! Return the zone index of the first no_mixing/anonymous_mixing zone
  ! scanning inward from center, limited to zones > lower_limit.
  ! The convective core occupies zones [result+1 .. nz].
  integer function conv_core_boundary_zone(id, lower_limit, ierr)
    integer, intent(in) :: id, lower_limit
    integer, intent(out) :: ierr
    type (star_info), pointer :: s
    integer :: zone
    ierr = 0
    call star_ptr(id, s, ierr)
    if (ierr /= 0) return
    do zone = s% nz, lower_limit + 1, -1
      if (s% mixing_type(zone) == no_mixing .or. &
          s% mixing_type(zone) == anonymous_mixing) then
        conv_core_boundary_zone = zone
        return
      end if
    end do
    conv_core_boundary_zone = s% nz
  end function conv_core_boundary_zone


  ! Return (as real) the zone index at the bottom of the outermost
  ! convective shell.  Returns real(nz+1) when no convective shell exists.
  real(dp) function conv_envelope_bottom_zone(id, ierr)
    integer, intent(in) :: id
    integer, intent(out) :: ierr
    type (star_info), pointer :: s
    integer :: zone, win_start, win_end, last_conv_zone, k_zone
    logical :: found_no_conv_window
    integer, parameter :: WINDOW = 100
    ierr = 0
    call star_ptr(id, s, ierr)
    if (ierr /= 0) return
    found_no_conv_window = .false.
    win_start = 1
    do while (win_start <= s% nz)
      win_end = min(win_start + WINDOW - 1, s% nz)
      if (.not. any(s% mixing_type(win_start:win_end) == convective_mixing)) then
        found_no_conv_window = .true.
        exit
      end if
      win_start = win_start + WINDOW
    end do
    k_zone = -1
    if (found_no_conv_window) then
      last_conv_zone = -1
      do zone = win_start - 1, 1, -1
        if (s% mixing_type(zone) == convective_mixing) then
          last_conv_zone = zone
          exit
        end if
      end do
      if (last_conv_zone > 0) then
        do zone = last_conv_zone + 1, s% nz
          if (s% mixing_type(zone) == no_mixing .or. &
              s% mixing_type(zone) == anonymous_mixing) then
            k_zone = zone
            exit
          end if
        end do
      end if
    end if
    if (k_zone == -1) then
      conv_envelope_bottom_zone = real(s% nz + 1, dp)
      return
    end if
    conv_envelope_bottom_zone = real(k_zone, dp)
  end function conv_envelope_bottom_zone


  subroutine data_for_extra_history_columns(id, id_extra, n, names, vals, ierr)
    integer, intent(in) :: id, id_extra, n
    character (len=maxlen_history_column_name) :: names(n)
    real(dp) :: vals(n)
    integer, intent(out) :: ierr
    type (star_info), pointer :: s

    integer :: i, k_core, conv_bottom_zone, env_bottom_zone
    real(dp) :: inorm, raw_env

    ierr = 0
    call star_ptr(id, s, ierr)
    if (ierr /= 0) return

    raw_env = conv_envelope_bottom_zone(id, ierr)
    if (int(raw_env) > s% nz) then
      env_bottom_zone = 1
    else
      env_bottom_zone = int(raw_env)
    end if

    k_core = conv_core_boundary_zone(id, env_bottom_zone, ierr)
    if (ierr /= 0) return
    conv_bottom_zone = k_core + 1

    if (conv_bottom_zone <= s% nz) then
      if (env_bottom_zone >= conv_bottom_zone) then
        env_bottom_zone = conv_bottom_zone
      end if
    end if

    inorm = 3.0d0 * Msun * Rsun * Rsun / 2.0d0

    names(1) = 'rcore'
    names(2) = 'mcore'
    if (conv_bottom_zone > s% nz) then
      vals(1) = 0.0d0
      vals(2) = 0.0d0
    else
      vals(1) = s% r(conv_bottom_zone) / Rsun
      vals(2) = s% m(conv_bottom_zone) / Msun
    end if

    names(3) = 'core_inertia'
    vals(3) = 0.0d0
    do i = conv_bottom_zone, s% nz
      vals(3) = vals(3) + s% dm(i) * s% rmid(i) * s% rmid(i)
    end do
    vals(3) = vals(3) / inorm

    names(4) = 'rrad'
    names(5) = 'mrad'
    if (env_bottom_zone <= 1 .or. env_bottom_zone > s% nz) then
      vals(4) = 0.0d0
      vals(5) = 0.0d0
    else
      vals(4) = s% r(env_bottom_zone) / Rsun
      vals(5) = s% m(env_bottom_zone) / Msun
    end if

    names(6) = 'rad_inertia'
    vals(6) = 0.0d0
    if (env_bottom_zone >= 1 .and. env_bottom_zone <= s% nz) then
      do i = env_bottom_zone, conv_bottom_zone - 1
        vals(6) = vals(6) + s% dm(i) * s% rmid(i) * s% rmid(i)
      end do
      vals(6) = vals(6) / inorm
    end if

    names(7) = 'env_inertia'
    vals(7) = 0.0d0
    do i = 1, env_bottom_zone - 1
      vals(7) = vals(7) + s% dm(i) * s% rmid(i) * s% rmid(i)
    end do
    vals(7) = vals(7) / inorm

    names(8) = 'core_outer_i'
    vals(8) = real(conv_bottom_zone, dp)

    names(9) = 'env_bottom_i'
    vals(9) = real(env_bottom_zone, dp)

  end subroutine data_for_extra_history_columns


  integer function how_many_extra_profile_columns(id, id_extra)
    use star_def, only: star_info
    integer, intent(in) :: id, id_extra
    integer :: ierr
    type (star_info), pointer :: s
    ierr = 0
    call star_ptr(id, s, ierr)
    if (ierr /= 0) return
    how_many_extra_profile_columns = 0
  end function how_many_extra_profile_columns


  subroutine data_for_extra_profile_columns(id, id_extra, n, nz, names, vals, ierr)
    use star_def, only: star_info, maxlen_profile_column_name
    integer, intent(in) :: id, id_extra, n, nz
    character (len=maxlen_profile_column_name) :: names(n)
    real(dp) :: vals(nz,n)
    integer, intent(out) :: ierr
    type (star_info), pointer :: s
    ierr = 0
    call star_ptr(id, s, ierr)
    if (ierr /= 0) return
  end subroutine data_for_extra_profile_columns


  integer function extras_finish_step(id, id_extra)
    integer, intent(in) :: id, id_extra
    integer :: ierr
    real(dp) :: envelope_mass_fraction, L_He, L_tot, min_center_h1_for_diff, &
         critmass, feh, rot_full_off, rot_full_on, frac2
    real(dp), parameter :: huge_dt_limit = 3.15d16
    real(dp), parameter :: new_varcontrol_target = 1d-3
    real(dp), parameter :: Zsol = 0.0142
    type (star_info), pointer :: s
    logical :: diff_test1, diff_test2, diff_test3
    character (len=strlen) :: photoname

    ierr = 0
    call star_ptr(id, s, ierr)
    if (ierr /= 0) return
    extras_finish_step = keep_going
    call store_extra_info(s)

    rot_full_off = s% job% extras_rpar(1)
    rot_full_on  = s% job% extras_rpar(2)

    if (rot_set_check) then
       if ((s% job% extras_rpar(3) > 0.0d0) .and. (s% initial_mass > rot_full_off)) then
          if ((abs(log10_cr(s% power_h_burn * Lsun / s% L(1))) < 1.0d-2) .and. (s% star_age > 1.0d2)) then
             if (s% initial_mass <= rot_full_on) then
                frac2 = (s% initial_mass - rot_full_off) / (rot_full_on - rot_full_off)
                frac2 = 0.5d0*(1.0d0 - cos(pi*frac2))
             else
                frac2 = 1.0d0
             end if
             s% job% new_omega_div_omega_crit = s% job% extras_rpar(3) * frac2
             s% job% set_near_zams_omega_div_omega_crit_steps = 10
             rot_set_check = .false.
          end if
       else
          rot_set_check = .false.
       end if
    end if

    if(TP_AGB_check .and. s% have_done_TP)then
       TP_AGB_check = .false.
       late_AGB_check = .true.
       call star_write_model(id, 'TPAGB.mod', ierr)
       photoname = 'photos/TPAGB_photo'
       call star_save_for_restart(id, photoname, ierr)
       s% delta_lgTeff_limit      = -1
       s% delta_lgTeff_hard_limit = -1
       s% delta_lgL_limit         = -1
       s% delta_lgL_hard_limit    = -1
       s% varcontrol_target = 2.0d0*s% varcontrol_target
    endif

    if(late_AGB_check)then
       if (s% initial_mass < 10.0d0 .and. s% he_core_mass/s% star_mass > 0.95d0) then
          photoname = 'photos/late_AGB_photo'
          call star_save_for_restart(id, photoname, ierr)
          late_AGB_check = .false.
          post_AGB_check = .true.
       endif
    endif

    if(post_AGB_check)then
       if(s% Teff > 1.0d4)then
          call star_write_model(id, 'postAGB.mod', ierr)
          photoname = 'photos/pAGB_photo'
          call star_save_for_restart(id, photoname, ierr)
          if(s% job% extras_lpar(2))then
             s% max_abar_for_burning = -1
          endif
          post_AGB_check = .false.
          pre_WD_check = .true.
          termination_code_str(t_xtra2) = 'stopping at post AGB'
          s% termination_code = t_xtra2
          extras_finish_step = terminate
       endif
    endif

    if(pre_WD_check)then
       if(s% Teff < 3.0d4 .and. s% L_surf < 1.0d0)then
          pre_WD_check = .false.
       endif
    endif

    ! Do not stop at TAMS here.  The MESA max_age control in the inlist
    ! determines the requested endpoint; the grid needs post-TAMS rows for
    ! systems whose final age is later than core-H depletion.

    if ((s% center_h1 < 1d-4) .and. (s% center_he4 < 5.0d-2) .and. (s% center_c12 < 5.0d-2)) then
       termination_code_str(t_xtra2) = 'central C12 mass fraction below 5e-2'
       s% termination_code = t_xtra2
       extras_finish_step = terminate
    endif

    if ((s% center_h1 < 1d-2) .and. (s% initial_mass <= 0.6d0) .and. (s% star_age > 5.0d10)) then
       termination_code_str(t_xtra2) = 'central H1 mass fraction below 0.01'
       s% termination_code = t_xtra2
       extras_finish_step = terminate
    endif

    min_center_h1_for_diff = 1d-10
    diff_test1 = abs(s% mass_conv_core - s% star_mass) < 1d-2
    diff_test2 = s% star_age > 1.0d11
    diff_test3 = s% center_h1 < min_center_h1_for_diff
    if (diff_test1 .or. diff_test2 .or. diff_test3) then
       s% diffusion_dt_limit = huge_dt_limit
    else
       s% diffusion_dt_limit = original_diffusion_dt_limit
    end if

  end function extras_finish_step


  subroutine Reimers_then_VW(id, Lsurf, Msurf, Rsurf, Tsurf, w, ierr)
    use star_def
    use chem_def, only: ih1, ihe4
    integer, intent(in) :: id
    real(dp), intent(in) :: Lsurf, Msurf, Rsurf, Tsurf
    real(dp), intent(out) :: w
    integer, intent(out) :: ierr
    real(dp) :: logP, P, reimers_w, pre_superwind_w, superwind_w, vexp, agb_w
    real(dp) :: center_h1, center_he4
    integer :: h1, he4
    type (star_info), pointer :: s
    ierr = 0
    call star_ptr(id, s, ierr)
    if (ierr /= 0) return
    reimers_w = 4d-13*(Lsurf*Rsurf/Msurf)/(Lsun*Rsun/Msun) * s% Reimers_scaling_factor
    logP = -2.07 + 1.94*log10_cr(Rsurf/Rsun) - 0.9*log10_cr(Msurf/Msun)
    P = pow_cr(10d0, logP)
    if (Msurf/Msun < 2.5) then
       pre_superwind_w = pow_cr(10d0, -11.4+0.0123*P)
    else
       pre_superwind_w = pow_cr(10d0, -11.4+0.0125*(P - 100.0*(Msurf/Msun - 2.5)))
    end if
    vexp = min(15.0, max(3.0, -13.5+0.056*P))*1e5
    superwind_w = (Lsurf/(clight*vexp))*(secyer/Msun)
    agb_w = min(pre_superwind_w, superwind_w)
    h1 = s% net_iso(ih1)
    he4 = s% net_iso(ihe4)
    center_h1  = s% xa(h1, s% nz)
    center_he4 = s% xa(he4, s% nz)
    if (center_h1 < 0.01d0 .and. center_he4 < s% RGB_to_AGB_wind_switch) then
       w = agb_w
    else
       w = reimers_w
    end if
  end subroutine Reimers_then_VW


  subroutine Reimers_then_Blocker(id, Lsurf, Msurf, Rsurf, Tsurf, w, ierr)
    use star_def
    use chem_def, only: ih1, ihe4
    integer, intent(in) :: id
    real(dp), intent(in) :: Lsurf, Msurf, Rsurf, Tsurf
    real(dp), intent(out) :: w
    integer, intent(out) :: ierr
    integer :: h1, he4
    real(dp) :: plain_reimers, reimers_w, blocker_w, center_h1, center_he4
    type (star_info), pointer :: s
    ierr = 0
    call star_ptr(id, s, ierr)
    if (ierr /= 0) return
    plain_reimers = 4d-13*(Lsurf*Rsurf/Msurf)/(Lsun*Rsun/Msun)
    reimers_w = plain_reimers * s% Reimers_scaling_factor
    blocker_w = plain_reimers * s% Blocker_scaling_factor * &
                4.83d-9 * pow_cr(Msurf/Msun,-2.1d0) * pow_cr(Lsurf/Lsun,2.7d0)
    h1  = s% net_iso(ih1)
    he4 = s% net_iso(ihe4)
    center_h1  = s% xa(h1, s% nz)
    center_he4 = s% xa(he4, s% nz)
    if (center_h1 < 0.01d0 .and. center_he4 > 0.1d0) then
       w = reimers_w
    else
       w = max(reimers_w, blocker_w)
    end if
  end subroutine Reimers_then_Blocker


  subroutine extras_after_evolve(id, id_extra, ierr)
    integer, intent(in) :: id, id_extra
    integer, intent(out) :: ierr
    type (star_info), pointer :: s
    ierr = 0
    call star_ptr(id, s, ierr)
    if (ierr /= 0) return
  end subroutine extras_after_evolve


  subroutine alloc_extra_info(s)
    integer, parameter :: extra_info_alloc = 1
    type (star_info), pointer :: s
    call move_extra_info(s, extra_info_alloc)
  end subroutine alloc_extra_info

  subroutine unpack_extra_info(s)
    integer, parameter :: extra_info_get = 2
    type (star_info), pointer :: s
    call move_extra_info(s, extra_info_get)
  end subroutine unpack_extra_info

  subroutine store_extra_info(s)
    integer, parameter :: extra_info_put = 3
    type (star_info), pointer :: s
    call move_extra_info(s, extra_info_put)
  end subroutine store_extra_info

  subroutine move_extra_info(s, op)
    integer, parameter :: extra_info_alloc = 1
    integer, parameter :: extra_info_get = 2
    integer, parameter :: extra_info_put = 3
    type (star_info), pointer :: s
    integer, intent(in) :: op
    integer :: i, num_ints, num_dbls, ierr
    i = 0
    num_ints = i
    i = 0
    num_dbls = i
    if (op /= extra_info_alloc) return
    if (num_ints == 0 .and. num_dbls == 0) return
    ierr = 0
    call star_alloc_extras(s% id, num_ints, num_dbls, ierr)
    if (ierr /= 0) then
       write(*,*) 'failed in star_alloc_extras'
       stop 1
    end if
  contains
    subroutine move_dbl(dbl)
      real(dp) :: dbl
      i = i+1
      select case (op)
      case (extra_info_get); dbl = s% extra_work(i)
      case (extra_info_put); s% extra_work(i) = dbl
      end select
    end subroutine move_dbl
    subroutine move_int(int)
      integer :: int
      i = i+1
      select case (op)
      case (extra_info_get); int = s% extra_iwork(i)
      case (extra_info_put); s% extra_iwork(i) = int
      end select
    end subroutine move_int
    subroutine move_flg(flg)
      logical :: flg
      i = i+1
      select case (op)
      case (extra_info_get); flg = (s% extra_iwork(i) /= 0)
      case (extra_info_put)
        if (flg) then
          s% extra_iwork(i) = 1
        else
          s% extra_iwork(i) = 0
        end if
      end select
    end subroutine move_flg
  end subroutine move_extra_info


  subroutine TSF(id, ierr)
    integer, intent(in) :: id
    integer, intent(out) :: ierr
    type (star_info), pointer :: s
    integer :: k, j, op_err, nsmooth, nsmootham
    real(dp) :: alpha, shearsmooth, nu_tsf, nu_tsf_t, omegac, omegag, omegaa, &
         omegat, difft, diffm, brunts, bruntsn2, logamnuomega, alphaq, scale_factor

    call star_ptr(id, s, ierr)
    if (ierr /= 0) return

    if(s% initial_mass > 4.0d0) then
       scale_factor = pow2((7d0 - s% initial_mass)/(7d0 - 4d0))
    else
       scale_factor = 1d0
    end if

    alpha = 1d0
    nsmooth = 5
    nsmootham = nsmooth-3
    shearsmooth = 1d-30
    op_err = 0

    do k = nsmooth+1, s% nz-(nsmooth+1)
       nu_tsf   = 1d-30
       nu_tsf_t = 1d-30
       shearsmooth = s% omega_shear(k)/(2d0*nsmooth+1d0)
       do j = 1, nsmooth
          shearsmooth = shearsmooth + (1d0/(2d0*nsmooth+1d0))* &
               (s% omega_shear(k-j) + s% omega_shear(k+j))
       end do
       diffm = diffmag(s% rho(k), s% T(k), s% abar(k), s% zbar(k), op_err)
       difft = 16d0*5.67d-5*pow3(s% T(k))/ &
               (3d0*s% opacity(k)*pow2(s% rho(k))*s% Cv(k))
       omegaa = s% omega(k)*pow_cr(shearsmooth*s% omega(k)/sqrt_cr(abs(s% brunt_N2(k))),1d0/3d0)
       omegat = difft*pow2(sqrt_cr(abs(s% brunt_N2(k)))/(omegaa*s% r(k)))
       brunts  = sqrt_cr(abs(s% brunt_N2_composition_term(k) + &
                   (s% brunt_N2(k)-s% brunt_N2_composition_term(k))/(1d0+omegat/omegaa)))
       bruntsn2 = sqrt_cr(abs(s% brunt_N2_composition_term(k) + &
                   (s% brunt_N2(k)-s% brunt_N2_composition_term(k))*min(1d0,diffm/difft)))
       brunts = max(brunts, bruntsn2)
       brunts = max(s% omega(k), brunts)
       omegaa = s% omega(k)*pow_cr(abs(shearsmooth*s% omega(k)/brunts),1d0/3d0)
       if (s% brunt_N2(k) > 0d0) then
          if (pow2(brunts) > 2d0*pow2(shearsmooth)*pow2(s% omega(k))) then
             omegac = 1d0*s% omega(k)*sqrt_cr(brunts/s% omega(k))* &
                      pow_cr(diffm/(pow2(s% r(k))*s% omega(k)),0.25d0)
             nu_tsf = 5d-1+5d-1*tanh(5d0*log(alpha*omegaa/omegac))
             nu_tsf = nu_tsf*pow3(alpha)*s% omega(k)*pow2(s% r(k))*pow2(s% omega(k)/brunts)
          end if
          if (pow2(brunts) < 2d0*pow2(shearsmooth)*pow2(s% omega(k))) then
             nu_tsf_t = alpha*abs(shearsmooth)*s% omega(k)*pow2(s% r(k))
          end if
          s% am_nu_omega(k) = s% am_nu_omega(k) + scale_factor*(max(nu_tsf,nu_tsf_t) + 1d-1)
       end if
    end do

    do k = s% nz-nsmooth, s% nz
       nu_tsf   = 1d-30
       nu_tsf_t = 1d-30
       diffm = diffmag(s% rho(k), s% T(k), s% abar(k), s% zbar(k), op_err)
       difft = 16d0*5.67d-5*pow3(s% T(k))/ &
               (3d0*s% opacity(k)*pow2(s% rho(k))*s% Cv(k))
       omegaa = s% omega(k)*pow_cr(shearsmooth*s% omega(k)/sqrt_cr(abs(s% brunt_N2(k))),1d0/3d0)
       omegat = difft*pow2(sqrt_cr(abs(s% brunt_N2(k)))/(s% omega(k)*s% r(k)))
       brunts  = sqrt_cr(abs(s% brunt_N2_composition_term(k) + &
                   (s% brunt_N2(k)-s% brunt_N2_composition_term(k))/(1d0+omegat/omegaa)))
       bruntsn2 = sqrt_cr(abs(s% brunt_N2_composition_term(k) + &
                   (s% brunt_N2(k)-s% brunt_N2_composition_term(k))*min(1d0,diffm/difft)))
       brunts = max(brunts, bruntsn2)
       brunts = max(s% omega(k), brunts)
       omegaa = s% omega(k)*pow_cr(abs(shearsmooth*s% omega(k)/brunts),1d0/3d0)
       if (s% brunt_N2(k) > 0d0) then
          if (pow2(brunts) > 2d0*pow2(shearsmooth)*pow2(s% omega(k))) then
             omegac = 1d0*s% omega(k)*sqrt_cr(brunts/s% omega(k))* &
                      pow_cr(diffm/(pow2(s% r(k))*s% omega(k)),0.25d0)
             nu_tsf = 5d-1+5d-1*tanh(5d0*log(alpha*omegaa/omegac))
             nu_tsf = nu_tsf*pow3(alpha)*s% omega(k)*pow2(s% r(k))*pow2(s% omega(k)/brunts)
          end if
          if (pow2(brunts) < 2d0*pow2(shearsmooth)*pow2(s% omega(k))) then
             nu_tsf_t = alpha*abs(shearsmooth)*s% omega(k)*pow2(s% r(k))
          end if
          s% am_nu_omega(k) = s% am_nu_omega(k) + scale_factor*(max(nu_tsf,nu_tsf_t) + 1d-1)
       end if
    end do

    do k = nsmooth, 1, -1
       nu_tsf   = 1d-30
       nu_tsf_t = 1d-30
       diffm = diffmag(s% rho(k), s% T(k), s% abar(k), s% zbar(k), op_err)
       difft = 16d0*5.67d-5*pow3(s% T(k))/ &
               (3d0*s% opacity(k)*pow2(s% rho(k))*s% Cv(k))
       omegaa = s% omega(k)*pow_cr(shearsmooth*s% omega(k)/sqrt_cr(abs(s% brunt_N2(k))),1d0/3d0)
       omegat = difft*pow2(sqrt_cr(abs(s% brunt_N2(k)))/(s% omega(k)*s% r(k)))
       brunts  = sqrt_cr(abs(s% brunt_N2_composition_term(k) + &
                   (s% brunt_N2(k)-s% brunt_N2_composition_term(k))/(1d0+omegat/omegaa)))
       bruntsn2 = sqrt_cr(abs(s% brunt_N2_composition_term(k) + &
                   (s% brunt_N2(k)-s% brunt_N2_composition_term(k))*min(1d0,diffm/difft)))
       brunts = max(brunts, bruntsn2)
       brunts = max(s% omega(k), brunts)
       omegaa = s% omega(k)*pow_cr(abs(shearsmooth*s% omega(k)/brunts),1d0/3d0)
       if (s% brunt_N2(k) > 0d0) then
          if (pow2(brunts) > 2d0*pow2(shearsmooth)*pow2(s% omega(k))) then
             omegac = 1d0*s% omega(k)*sqrt_cr(brunts/s% omega(k))* &
                      pow_cr(diffm/(pow2(s% r(k))*s% omega(k)),0.25d0)
             nu_tsf = 5d-1+5d-1*tanh(5d0*log(alpha*omegaa/omegac))
             nu_tsf = nu_tsf*pow3(alpha)*s% omega(k)*pow2(s% r(k))*pow2(s% omega(k)/brunts)
          end if
          if (pow2(brunts) < 2d0*pow2(shearsmooth)*pow2(s% omega(k))) then
             nu_tsf_t = alpha*abs(shearsmooth)*s% omega(k)*pow2(s% r(k))
          end if
          s% am_nu_omega(k) = s% am_nu_omega(k) + scale_factor*(max(nu_tsf,nu_tsf_t) + 1d-1)
       end if
    end do

    logamnuomega = -3d1
    do k = nsmootham+1, s% nz-(nsmootham+1)
       if (s% mixing_type(k) == 1) then
          s% am_nu_omega(k) = s% am_nu_omega(k)
       else
          logamnuomega = log10_cr(s% am_nu_omega(k))/(2d0*nsmootham+1d0)
       end if
       do j = 1, nsmootham
          if (s% mixing_type(k-j) < 3.5d0) then
             logamnuomega = log10_cr(s% am_nu_omega(k))
          else
             logamnuomega = logamnuomega + (1d0/(2d0*nsmootham+1d0))*log10_cr(s% am_nu_omega(k-j))
          end if
       end do
       do j = 1, nsmootham
          if (s% mixing_type(k+j) < 3.5d0) then
             logamnuomega = logamnuomega
          else
             logamnuomega = logamnuomega + (1d0/(2d0*nsmootham+1d0))*log10_cr(s% am_nu_omega(k+j))
          end if
       end do
       s% am_nu_omega(k) = exp10_cr(logamnuomega)
    end do

    do k = s% nz-nsmootham, s% nz
       s% am_nu_omega(k) = s% am_nu_omega(k-1)
    end do
    do k = nsmootham, 1, -1
       s% am_nu_omega(k) = s% am_nu_omega(k-1)
    end do

  end subroutine TSF


  real(dp) function diffmag(rho, T, abar, zbar, ierr)
    real(dp), intent(in) :: rho, T, abar, zbar
    integer, intent(out) :: ierr
    real(dp) :: xgamma, xlg, xsig1, xsig2, xsig3, xxx, ffff, xsig

    if (ierr /= 0) return

    xgamma = 0.2275d0*zbar*zbar*pow_cr(rho*1d-6/abar,1d0/3d0)*1d8/T
    xlg = log10_cr(xgamma)
    if (xlg < -1.5d0) then
       xsig = sige1(zbar, T, xgamma)
    else if (xlg >= -1.5d0 .and. xlg <= 0d0) then
       xxx  = (xlg + 0.75d0)*4d0/3d0
       ffff = 0.25d0*(2d0 - 3d0*xxx + xxx*xxx*xxx)
       xsig1 = sige1(zbar, T, xgamma)
       xsig2 = sige2(T, rho, zbar, ierr)
       if (ierr /= 0) return
       xsig = (1d0-ffff)*xsig2 + ffff*xsig1
    else if (xlg > 0d0 .and. xlg < 0.5d0) then
       xsig = sige2(T, rho, zbar, ierr)
       if (ierr /= 0) return
    else if (xlg >= 0.5d0 .and. xlg < 1d0) then
       xxx  = (xlg - 0.75d0)*4d0
       ffff = 0.25d0*(2d0 - 3d0*xxx + xxx*xxx*xxx)
       xsig2 = sige2(T, rho, zbar, ierr)
       if (ierr /= 0) return
       xsig3 = sige3(zbar, T, xgamma)
       xsig = (1d0-ffff)*xsig3 + ffff*xsig2
    else
       xsig = sige3(zbar, T, xgamma)
    endif
    diffmag = 7.1520663d19/xsig

  end function diffmag


  real(dp) function sqrt_cr(x)
    real(dp), intent(in) :: x
    sqrt_cr = pow_cr(x, 0.5d0)
  end function sqrt_cr


  real(dp) function sige1(z, t, xgamma)
    real(dp), intent(in) :: z, t, xgamma
    real(dp) :: etan, xlambda, f
    if (t >= 4.2d5) then
       f = sqrt_cr(4.2d5/t)
    else
       f = 1d0
    end if
    xlambda = sqrt_cr(3d0*z*z*z)*pow_cr(xgamma,-1.5d0)*f + 1d0
    etan = 3d11*z*log_cr(xlambda)*pow_cr(t,-1.5d0)
    etan = etan/(1d0-1.20487d0*exp_cr(-1.0576d0*pow_cr(z,0.347044d0)))
    sige1 = clight*clight/(4d0*pi*etan)
  end function sige1


  real(dp) function sige2(T, rho, zbar, ierr)
    use kap_lib, only: kap_get_elect_cond_opacity
    real(dp), intent(in) :: T, rho, zbar
    integer, intent(out) :: ierr
    real(dp) :: kap, dlnkap_dlnRho, dlnkap_dlnT
    call kap_get_elect_cond_opacity( &
         zbar, log10_cr(rho), log10_cr(T), &
         kap, dlnkap_dlnRho, dlnkap_dlnT, ierr)
    sige2 = 1.11d9*T*T/(rho*kap)
  end function sige2


  real(dp) function sige3(z, t, xgamma)
    real(dp), intent(in) :: z, t, xgamma
    real(dp) :: rme, rm23, ctmp, xi
    rme = 8.5646d-23*t*t*t*xgamma*xgamma*xgamma/powi_cr(z,5)
    rm23 = pow_cr(rme, 2d0/3d0)
    ctmp = 1d0 + 1.018d0*rm23
    xi = sqrt_cr(3.14159d0/3d0)*log_cr(z)/3d0 + &
         2d0*log_cr(1.32d0+2.33d0/sqrt_cr(xgamma))/3d0 - 0.484d0*rm23/ctmp
    sige3 = 8.630d21*rme/(z*ctmp*xi)
  end function sige3


end module run_star_extras

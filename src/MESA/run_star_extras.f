! ***********************************************************************
!
!   Copyright (C) 2010  Bill Paxton
!
!   this file is part of mesa.
!
!   mesa is free software; you can redistribute it and/or modify
!   it under the terms of the gnu general library public license as published
!   by the free software foundation; either version 2 of the license, or
!   (at your option) any later version.
!
!   mesa is distributed in the hope that it will be useful, 
!   but without any warranty; without even the implied warranty of
!   merchantability or fitness for a particular purpose.  see the
!   gnu library general public license for more details.
!
!   you should have received a copy of the gnu library general public license
!   along with this software; if not, write to the free software
!   foundation, inc., 59 temple place, suite 330, boston, ma 02111-1307 usa
!
! ***********************************************************************
 
      module run_star_extras

      use star_lib
      use star_def
      use const_def
      
      implicit none
      
      ! these routines are called by the standard run_star check_model
      contains
      
      subroutine extras_controls(s, ierr)
         type (star_info), pointer :: s
         integer, intent(out) :: ierr
         ierr = 0
         
         ! this is the place to set any procedure pointers you want to change
         ! e.g., other_wind, other_mixing, other_energy  (see star_data.inc)
         
         
      end subroutine extras_controls
      
      
      integer function extras_startup(s, id, restart, ierr)
         type (star_info), pointer :: s
         integer, intent(in) :: id
         logical, intent(in) :: restart
         integer, intent(out) :: ierr
         ierr = 0
         extras_startup = 0
         if (.not. restart) then
            call alloc_extra_info(s)
         else ! it is a restart
            call unpack_extra_info(s)
         end if
      end function extras_startup
      

      ! returns either keep_going, retry, backup, or terminate.
      integer function extras_check_model(s, id, id_extra)
         type (star_info), pointer :: s
         integer, intent(in) :: id, id_extra
         extras_check_model = keep_going         
         if (.false. .and. s% star_mass_h1 < 0.35d0) then
            ! stop when star hydrogen mass drops to specified level
            extras_check_model = terminate
            write(*, *) 'have reached desired hydrogen mass'
            return
         end if
      end function extras_check_model


      integer function how_many_extra_history_columns(s, id, id_extra)
         type (star_info), pointer :: s
         integer, intent(in) :: id, id_extra
         !2013 Oct 14 -- Adding Irot column
         how_many_extra_history_columns = 2
      end function how_many_extra_history_columns
      
      
      subroutine data_for_extra_history_columns(s, id, id_extra, n, names, vals, ierr)
         type (star_info), pointer :: s
         integer, intent(in) :: id, id_extra, n
         character (len=maxlen_history_column_name) :: names(n)
         real(dp) :: vals(n)
         !cgs units
         real(dp), parameter :: MSun = 1.99d33, RSun = 6.955d10
         integer, intent(out) :: ierr

         integer :: i
         !moments of inertia for radiative and convective zones
         real(dp) :: bot_r, top_r, &
                     conv_mx1_mid_r, conv_mx2_mid_r, &
                     conv_mx_bot_r, conv_mx_top_r, &
                     Irot_rad, Irot_conv
         
         !note: do NOT add these names to history_columns.list
         ! the history_columns.list is only for the built-in log column options.
         ! it must not include the new column names you are adding here.

         !2013 Nov 6 -- N.B. This code ASSUMES that the convective zone lies
         !  exterior to the radiative zone. Make sure that assumption works
         !  before using this code.
         Irot_conv = 0.
         Irot_rad = 0.
     
         !s%r(i) is the radius at the OUTER edge of cell i
         !
         !Just to confuse things, the bottom cell is the i = nz cell, 
         !  while the top cell is the i = 1 cell.
         do i = s%nz - 1, 1, -1 

           !the top and bottom of the current cell -- 

           !the top and bottom of the current cell -- 
           !  Note that bottomost (nz-th) cell has its bottom at r = 0.
           bot_r = s%r(i+1)/RSun
           top_r = s%r(i)/RSun

           !2014 Apr 11 -- Consider the outermost CZ, not the biggest
           conv_mx1_mid_r = 0.5*(s%conv_mx1_top_r + s%conv_mx1_bot_r)
           conv_mx2_mid_r = 0.5*(s%conv_mx2_top_r + s%conv_mx2_bot_r)

           conv_mx_top_r = s%conv_mx1_top_r
           conv_mx_bot_r = s%conv_mx1_bot_r

           if(conv_mx2_mid_r > conv_mx1_mid_r) then 
             conv_mx_top_r = s%conv_mx2_top_r
             conv_mx_bot_r = s%conv_mx2_bot_r
           end if

           !If bottom of the current cell is above the bottom of the conv AND
           !                                 below the top                OR
           !   top    of the current cell is below the top                AND
           !                                 above the bottom
           if(( (bot_r >= conv_mx_bot_r) .and. &
                (bot_r <= conv_mx_top_r) ) & 
              .or. &
              ( (top_r <= conv_mx_top_r) .and. &
                (top_r >= conv_mx_bot_r) ) ) then

              !shift boundaries over which to integrate, if necessary
              if(conv_mx_bot_r > bot_r) bot_r = conv_mx_bot_r
              if(conv_mx_top_r < top_r) top_r = conv_mx_top_r

              bot_r = bot_r*Rsun
              top_r = top_r*Rsun

              Irot_conv = Irot_conv + 2./5*s%dm(i)* &
                (bot_r**4 + bot_r*top_r*&
                  (bot_r**2 + bot_r*top_r + top_r**2) + &
                 top_r**4)/&
                (bot_r**2 + bot_r*top_r + top_r**2)
           end if

           if(bot_r < conv_mx_bot_r) then 

              !shift boundaries over which to integrate, if necessary
              if(top_r > conv_mx_bot_r) top_r = conv_mx_bot_r

              bot_r = bot_r*Rsun
              top_r = top_r*Rsun

              Irot_rad = Irot_rad + 2./5*s%dm(i)* &
                (bot_r**4 + bot_r*top_r*&
                  (bot_r**2 + bot_r*top_r + top_r**2) + &
                 top_r**4)/&
                (bot_r**2 + bot_r*top_r + top_r**2)
           end if
        
         end do

         names(1) = "Irot_conv"
         vals(1) = Irot_conv
         
         names(2) = "Irot_rad"
         vals(2) = Irot_rad

         ierr = 0
      end subroutine data_for_extra_history_columns

      
      integer function how_many_extra_profile_columns(s, id, id_extra)
         type (star_info), pointer :: s
         integer, intent(in) :: id, id_extra
         how_many_extra_profile_columns = 0
      end function how_many_extra_profile_columns
      
      
      subroutine data_for_extra_profile_columns(s, id, id_extra, n, nz, names, vals, ierr)
         type (star_info), pointer :: s
         integer, intent(in) :: id, id_extra, n, nz
         character (len=maxlen_profile_column_name) :: names(n)
         real(dp) :: vals(nz,n)
         integer, intent(out) :: ierr
         integer :: k
         ierr = 0
         
         !note: do NOT add these names to profile_columns.list
         ! the profile_columns.list is only for the built-in profile column options.
         ! it must not include the new column names you are adding here.

         ! here is an example for adding a profile column
         !if (n /= 1) stop 'data_for_extra_profile_columns'
         !names(1) = 'beta'
         !do k = 1, nz
         !   vals(k,1) = s% Pgas(k)/s% P(k)
         !end do
         
      end subroutine data_for_extra_profile_columns
      

      ! returns either keep_going or terminate.
      ! note: cannot request retry or backup; extras_check_model can do that.
      integer function extras_finish_step(s, id, id_extra)
         type (star_info), pointer :: s
         integer, intent(in) :: id, id_extra
         integer :: ierr
         extras_finish_step = keep_going
         call store_extra_info(s)

         ! to save a profile, 
            ! s% need_to_save_profiles_now = .true.
         ! to update the star log,
            ! s% need_to_update_history_now = .true.
            
      end function extras_finish_step
      
      
      subroutine extras_after_evolve(s, id, id_extra, ierr)
         type (star_info), pointer :: s
         integer, intent(in) :: id, id_extra
         integer, intent(out) :: ierr
         ierr = 0
      end subroutine extras_after_evolve
      
      
      ! routines for saving and restoring extra data so can do restarts
         
         ! put these defs at the top and delete from the following routines
         !integer, parameter :: extra_info_alloc = 1
         !integer, parameter :: extra_info_get = 2
         !integer, parameter :: extra_info_put = 3
      
      
      subroutine alloc_extra_info(s)
         integer, parameter :: extra_info_alloc = 1
         type (star_info), pointer :: s
         call move_extra_info(s,extra_info_alloc)
      end subroutine alloc_extra_info
      
      
      subroutine unpack_extra_info(s)
         integer, parameter :: extra_info_get = 2
         type (star_info), pointer :: s
         call move_extra_info(s,extra_info_get)
      end subroutine unpack_extra_info
      
      
      subroutine store_extra_info(s)
         integer, parameter :: extra_info_put = 3
         type (star_info), pointer :: s
         call move_extra_info(s,extra_info_put)
      end subroutine store_extra_info
      
      
      subroutine move_extra_info(s,op)
         integer, parameter :: extra_info_alloc = 1
         integer, parameter :: extra_info_get = 2
         integer, parameter :: extra_info_put = 3
         type (star_info), pointer :: s
         integer, intent(in) :: op
         
         integer :: i, j, num_ints, num_dbls, ierr
         
         i = 0
         ! call move_int or move_flg    
         num_ints = i
         
         i = 0
         ! call move_dbl       
         
         num_dbls = i
         
         if (op /= extra_info_alloc) return
         if (num_ints == 0 .and. num_dbls == 0) return
         
         ierr = 0
         call star_alloc_extras(s% id, num_ints, num_dbls, ierr)
         if (ierr /= 0) then
            write(*,*) 'failed in star_alloc_extras'
            write(*,*) 'alloc_extras num_ints', num_ints
            write(*,*) 'alloc_extras num_dbls', num_dbls
            stop 1
         end if
         
         contains
         
         subroutine move_dbl(dbl)
            real(dp) :: dbl
            i = i+1
            select case (op)
            case (extra_info_get)
               dbl = s% extra_work(i)
            case (extra_info_put)
               s% extra_work(i) = dbl
            end select
         end subroutine move_dbl
         
         subroutine move_int(int)
            integer :: int
            i = i+1
            select case (op)
            case (extra_info_get)
               int = s% extra_iwork(i)
            case (extra_info_put)
               s% extra_iwork(i) = int
            end select
         end subroutine move_int
         
         subroutine move_flg(flg)
            logical :: flg
            i = i+1
            select case (op)
            case (extra_info_get)
               flg = (s% extra_iwork(i) /= 0)
            case (extra_info_put)
               if (flg) then
                  s% extra_iwork(i) = 1
               else
                  s% extra_iwork(i) = 0
               end if
            end select
         end subroutine move_flg
      
      end subroutine move_extra_info
      end module run_star_extras
      

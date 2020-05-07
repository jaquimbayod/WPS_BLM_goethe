!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! MODULE SMOOTH_MODULE
!
! This module provides routines for smoothing.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module smooth_module

   use parallel_module
!B: Added by CC
   use gridinfo_module
   use misc_definitions_module
   use module_debug
!E: Added by CC

   contains

 
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! Name: one_two_one
   !
   ! Purpose: Apply the 1-2-1 smoother from the MM5 program TERRAIN 
   !   (found in smth121.F) to array.
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine one_two_one(array, start_dom_x, end_dom_x, start_dom_y, end_dom_y, &
                          start_x, end_x, start_y, end_y, start_z, end_z, npass, msgval)
 
      implicit none
  
      ! Arguments
      integer, intent(in) :: start_dom_x, start_dom_y, start_x, start_y, start_z
      integer, intent(in) :: end_dom_x, end_dom_y, end_x, end_y, end_z
      integer, intent(in) :: npass
      real, intent(in) :: msgval
      real, dimension(start_x:end_x, start_y:end_y, start_z:end_z), intent(inout) :: array 
  
      ! Local variables
      integer :: ix, iy, iz, ipass
      real, pointer, dimension(:,:,:) :: scratch
  
      allocate(scratch(start_x+1:end_x-1, start_y:end_y, start_z:end_z))
  
      do ipass=1,npass

         do iy=start_y,end_y
            do ix=start_x+1,end_x-1
               do iz=start_z,end_z
                  scratch(ix,iy,iz) = 0.50*array(ix,iy,iz)+0.25*(array(ix-1,iy,iz)+array(ix+1,iy,iz))
               end do
            end do
         end do
   
         do iy=start_y+1,end_y-1
            do ix=start_x+1,end_x-1
               do iz=start_z,end_z
                  array(ix,iy,iz) = 0.50*scratch(ix,iy,iz)+0.25*(scratch(ix,iy-1,iz)+scratch(ix,iy+1,iz))
               end do
             end do
          end do

         call exchange_halo_r(array, &
                              start_x, end_x, start_y, end_y, start_z, end_z, &
                              start_dom_x, end_dom_x, start_dom_y, end_dom_y, start_z, end_z)

      end do
  
      deallocate(scratch)
 
   end subroutine one_two_one 
 
 
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! Name: smth_desmth
   !
   ! Purpose: Apply the smoother-desmoother from the MM5 program TERRAIN 
   !   (found in smther.F) to array.
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine smth_desmth(array, start_dom_x, end_dom_x, start_dom_y, end_dom_y, &
                          start_x, end_x, start_y, end_y, start_z, end_z, npass, msgval)
 
      implicit none
  
      ! Arguments
      integer, intent(in) :: start_dom_x, start_dom_y, start_x, start_y, start_z
      integer, intent(in) :: end_dom_x, end_dom_y, end_x, end_y, end_z
      integer, intent(in) :: npass
      real, intent(in) :: msgval
      real, dimension(start_x:end_x, start_y:end_y, start_z:end_z), intent(inout) :: array 
  
      ! Local variables
      integer :: ix, iy, iz, ipass
      real, pointer, dimension(:,:,:) :: scratch
  
      allocate(scratch(start_x+1:end_x-1, start_y:end_y, start_z:end_z))
  
      do ipass=1,npass

         !
         ! Smoothing pass
         !
         do iy=start_y,end_y
            do ix=start_x+1,end_x-1
               do iz=start_z,end_z
                  scratch(ix,iy,iz) = 0.5*array(ix,iy,iz) + 0.25*(array(ix-1,iy,iz)+array(ix+1,iy,iz))
               end do
            end do
         end do
   
         do iy=start_y+1,end_y-1
            do ix=start_x+1,end_x-1
               do iz=start_z,end_z
                  array(ix,iy,iz) = 0.5*scratch(ix,iy,iz) + 0.25*(scratch(ix,iy-1,iz)+scratch(ix,iy+1,iz))
               end do
            end do
         end do

         call exchange_halo_r(array, &
                              start_x, end_x, start_y, end_y, start_z, end_z, &
                              start_dom_x, end_dom_x, start_dom_y, end_dom_y, start_z, end_z)
   
         !
         ! Desmoothing pass
         !
         do iy=start_y,end_y
            do ix=start_x+1,end_x-1
               do iz=start_z,end_z
                  scratch(ix,iy,iz) = 1.52*array(ix,iy,iz) - 0.26*(array(ix-1,iy,iz)+array(ix+1,iy,iz))
               end do
            end do
         end do
   
         do iy=start_y+1,end_y-1
            do ix=start_x+1,end_x-1
               do iz=start_z,end_z
                  array(ix,iy,iz) = 1.52*scratch(ix,iy,iz) - 0.26*(scratch(ix,iy-1,iz)+scratch(ix,iy+1,iz))
               end do
            end do
         end do

         call exchange_halo_r(array, &
                              start_x, end_x, start_y, end_y, start_z, end_z, &
                              start_dom_x, end_dom_x, start_dom_y, end_dom_y, start_z, end_z)

      end do
  
      deallocate(scratch)
 
   end subroutine smth_desmth


   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! Name: smth_desmth_special
   !
   ! Purpose: Apply the smoother-desmoother from the MM5 program TERRAIN 
   !   (found in smther.F) to array; however, any grid points that were not
   !   originally negative but which have been smoothed to a negative value
   !   will be restored to their original values.
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine smth_desmth_special(array, start_dom_x, end_dom_x, start_dom_y, end_dom_y, &
                                  start_x, end_x, start_y, end_y, start_z, end_z, npass, msgval)

      implicit none

      ! Arguments
      integer, intent(in) :: start_dom_x, start_dom_y, start_x, start_y, start_z
      integer, intent(in) :: end_dom_x, end_dom_y, end_x, end_y, end_z
      integer, intent(in) :: npass
      real, intent(in) :: msgval
      real, dimension(start_x:end_x, start_y:end_y, start_z:end_z), intent(inout) :: array

      ! Local variables
      integer :: ix, iy, iz, ipass
      real, pointer, dimension(:,:,:) :: scratch, orig_array

      allocate(scratch(start_x+1:end_x-1, start_y:end_y, start_z:end_z))
      allocate(orig_array(start_x:end_x, start_y:end_y, start_z:end_z))

      orig_array = array

      do ipass=1,npass

         !
         ! Smoothing pass
         !
         do iy=start_y,end_y
            do ix=start_x+1,end_x-1
               do iz=start_z,end_z
                  scratch(ix,iy,iz) = 0.5*array(ix,iy,iz) + 0.25*(array(ix-1,iy,iz)+array(ix+1,iy,iz))
               end do
            end do
         end do

         do iy=start_y+1,end_y-1
            do ix=start_x+1,end_x-1
               do iz=start_z,end_z
                  array(ix,iy,iz) = 0.5*scratch(ix,iy,iz) + 0.25*(scratch(ix,iy-1,iz)+scratch(ix,iy+1,iz))
               end do
            end do
         end do

         call exchange_halo_r(array, &
                              start_x, end_x, start_y, end_y, start_z, end_z, &
                              start_dom_x, end_dom_x, start_dom_y, end_dom_y, start_z, end_z)

         !
         ! Desmoothing pass
         !
         do iy=start_y,end_y
            do ix=start_x+1,end_x-1
               do iz=start_z,end_z
                  scratch(ix,iy,iz) = 1.52*array(ix,iy,iz) - 0.26*(array(ix-1,iy,iz)+array(ix+1,iy,iz))
               end do
            end do
         end do

         do iy=start_y+1,end_y-1
            do ix=start_x+1,end_x-1
               do iz=start_z,end_z
                  array(ix,iy,iz) = 1.52*scratch(ix,iy,iz) - 0.26*(scratch(ix,iy-1,iz)+scratch(ix,iy+1,iz))
               end do
            end do
         end do

         call exchange_halo_r(array, &
                              start_x, end_x, start_y, end_y, start_z, end_z, &
                              start_dom_x, end_dom_x, start_dom_y, end_dom_y, start_z, end_z)

      end do

      ! Remove artificially negative values
      do iy=start_y,end_y
         do ix=start_x,end_x
            do iz=start_z,end_z
               if (array(ix,iy,iz) < 0. .and. orig_array(ix,iy,iz) >= 0.) then
                  array(ix,iy,iz) = orig_array(ix,iy,iz)
               end if
            end do
         end do
      end do

      deallocate(scratch)
      deallocate(orig_array)

   end subroutine smth_desmth_special

!B: Added by CC
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! Name: smth_cc_ardu
   !
   ! Purpose: Reduce any angle above a specified threshold if required.
   !   Adapted from the subroutine smth_desmth_special.
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!  Pass which_domain
   subroutine smth_cc_ardu(which_domain, array, start_dom_x, end_dom_x, start_dom_y, end_dom_y, &
                           start_x, end_x, start_y, end_y, start_z, end_z, smth_opt, npass, msgval)

      implicit none

      ! Arguments
      integer, intent(in) :: start_dom_x, start_dom_y, start_x, start_y, start_z
      integer, intent(in) :: end_dom_x, end_dom_y, end_x, end_y, end_z
      integer, intent(in) :: npass
      real, intent(in) :: msgval
      real, dimension(start_x:end_x, start_y:end_y, start_z:end_z), intent(inout) :: array

      ! Local variables
      integer :: ix, iy, iz, ipass
      real, pointer, dimension(:,:,:) :: scratch, orig_array
      ! Add new variables
      integer, intent(in) :: which_domain
      integer, intent(in) :: smth_opt
      integer :: id
      real :: resolution_x, resolution_y, slope_angle_x_max, slope_angle_y_max, slope_angle_max
      real, allocatable, dimension(:,:,:) :: slope_angle_x, slope_angle_y
      logical, allocatable, dimension(:,:,:) :: do_it_x, do_it_x_u, do_it_x_d, do_it_y, do_it_y_u, do_it_y_d
      real, parameter :: epsilon_slope_angle = 0.1
      integer :: counter, number_pts
      integer, parameter :: max_counter = 100

      allocate(scratch(start_x:end_x, start_y:end_y, start_z:end_z))
      allocate(orig_array(start_x:end_x, start_y:end_y, start_z:end_z))
      ! Allocate new variables
      allocate(slope_angle_x(start_x:end_x, start_y:end_y, start_z:end_z))
      allocate(slope_angle_y(start_x:end_x, start_y:end_y, start_z:end_z))
      allocate(do_it_x(start_x:end_x, start_y:end_y, start_z:end_z))
      allocate(do_it_x_u(start_x:end_x, start_y:end_y, start_z:end_z))
      allocate(do_it_x_d(start_x:end_x, start_y:end_y, start_z:end_z))
      allocate(do_it_y(start_x:end_x, start_y:end_y, start_z:end_z))
      allocate(do_it_y_u(start_x:end_x, start_y:end_y, start_z:end_z))
      allocate(do_it_y_d(start_x:end_x, start_y:end_y, start_z:end_z))

      orig_array = array
!     Determine dx and dy for the domain
      resolution_x = dxkm
      resolution_y = dykm
      do id = 1, which_domain
        resolution_x = resolution_x / parent_grid_ratio(id)
        resolution_y = resolution_y / parent_grid_ratio(id)
      end do

!     Print target angle
      call mprintf(.true.,LOGFILE,'    Target angle: %f',f1=max_slope_angle(id-1))

!     Calculate slope_angle_max
!     In the x-direction
      do_it_x_u = .false.
      do_it_x_d = .false.
      number_pts = 0
      do iy=start_y+3,end_y-3
         do ix=start_x+3,end_x-3
            do iz=start_z,end_z
      	       slope_angle_x(ix,iy,iz) = abs(array(ix,iy,iz) - array(ix-1,iy,iz)) / resolution_x
      	       slope_angle_x(ix,iy,iz) = atan(slope_angle_x(ix,iy,iz)) * 180. / (4 * atan(1.))
      	       if (slope_angle_x(ix,iy,iz) > max_slope_angle(id-1) + epsilon_slope_angle) then 
      	     	  number_pts = number_pts + 1
      	     	  if (array(ix,iy,iz) > array(ix-1,iy,iz)) then
      	     	     do_it_x_u(ix,iy,iz) = .false.
      	     	     do_it_x_d(ix,iy,iz) = .true.
      	     	  else
      	     	     do_it_x_u(ix,iy,iz) = .true.
      	     	     do_it_x_d(ix,iy,iz) = .false.
      	     	  end if
      	       end if
            end do
         end do
      end do
      call mprintf(.true.,LOGFILE,'    %i point(s) exceed(s) the target angle (%f) in the x-direction',i1=number_pts,f1=max_slope_angle(id-1))
      if ((number_pts > 0) .and. (debug_level >= 100)) then
         do iy=start_y+3,end_y-3
            do ix=start_x+3,end_x-3
               do iz=start_z,end_z
        	  if (do_it_x_u(ix,iy,iz) .or. do_it_x_d(ix,iy,iz)) then
        	     call mprintf(.true.,LOGFILE,'	   at point: i=%i, j=%i',i1=ix+3,i2=iy+3)
        	  end if
               end do
            end do
         end do 	   
      end if
!     In the y-direction
      do_it_y_u = .false.
      do_it_y_d = .false.
      number_pts = 0
      do iy=start_y+3,end_y-3
         do ix=start_x+3,end_x-3
            do iz=start_z,end_z
      	       slope_angle_y(ix,iy,iz) = abs(array(ix,iy,iz) - array(ix,iy-1,iz)) / resolution_y
      	       slope_angle_y(ix,iy,iz) = atan(slope_angle_y(ix,iy,iz)) * 180. / (4 * atan(1.))
      	       if (slope_angle_y(ix,iy,iz) > max_slope_angle(id-1) + epsilon_slope_angle) then
      	     	  number_pts = number_pts + 1
      	     	  if (array(ix,iy,iz) > array(ix,iy-1,iz)) then
      	     	     do_it_y_u(ix,iy,iz) = .false.
      	     	     do_it_y_d(ix,iy,iz) = .true.
      	     	  else
      	     	     do_it_y_u(ix,iy,iz) = .true.
      	     	     do_it_y_d(ix,iy,iz) = .false.
      	     	  end if
      	       end if
            end do
         end do
      end do
      call mprintf(.true.,LOGFILE,'    %i point(s) exceed(s) the target angle (%f) in the y-direction',i1=number_pts,f1=max_slope_angle(id-1))
      if ((number_pts > 0) .and. (debug_level >= 100)) then
         do iy=start_y+3,end_y-3
            do ix=start_x+3,end_x-3
               do iz=start_z,end_z
        	  if (do_it_y_u(ix,iy,iz) .or. do_it_y_d(ix,iy,iz)) then
        	     call mprintf(.true.,LOGFILE,'	   at point: i=%i, j=%i',i1=ix+3,i2=iy+3)
        	  end if
               end do
            end do
         end do 	   
      end if
      number_pts = 0
      do iy=start_y+3,end_y-3
         do ix=start_x+3,end_x-3
            do iz=start_z,end_z
     	       if ((do_it_x_d(ix,iy,iz) .and. do_it_y_u(ix,iy,iz)) .or. &
     	     	   (do_it_x_u(ix,iy,iz) .and. do_it_y_d(ix,iy,iz))) then
     	     	  number_pts = number_pts + 1
     	       end if
            end do
         end do
      end do
      call mprintf(.true.,LOGFILE,'    %i saddle point(s) exceed(s) the target angle (%f) in both the x- and y-directions',i1=number_pts,f1=max_slope_angle(id-1))
      if ((number_pts > 0) .and. (debug_level >= 100)) then
         do iy=start_y+3,end_y-3
            do ix=start_x+3,end_x-3
               do iz=start_z,end_z
        	  if ((do_it_x_d(ix,iy,iz) .and. do_it_y_u(ix,iy,iz)) .or. &
        	      (do_it_x_u(ix,iy,iz) .and. do_it_y_d(ix,iy,iz))) then
        	     call mprintf(.true.,LOGFILE,'	   at point: i=%i, j=%i',i1=ix+3,i2=iy+3)
        	  end if
               end do
            end do
         end do 	   
      end if
      slope_angle_x_max = MAXVAL(slope_angle_x(start_x+3:end_x-3,start_y+3:end_y-3,:))
      slope_angle_y_max = MAXVAL(slope_angle_y(start_x+3:end_x-3,start_y+3:end_y-3,:))
      slope_angle_max = MAX(slope_angle_x_max,slope_angle_y_max)

      call mprintf(.true.,LOGFILE,'    Maximum angle: %f',f1=slope_angle_max)

!     Scale the topography if requested. Only once.

      if (scale_topo(id-1)) then

         call mprintf(.true.,LOGFILE,'    Scale the topography by %f as requested',f1=scale_topo_factor(id-1))

         do iy=start_y,end_y
            do ix=start_x,end_x
               do iz=start_z,end_z
        	  array(ix,iy,iz) = array(ix,iy,iz) * (1 - (1 - scale_topo_factor(id-1)) * & 
        						   (array(ix,iy,iz) - minval(array(:,:,iz))) / &
        						   (maxval(array(:,:,iz)) - minval(array(:,:,iz))))
               end do
            end do
         end do

!        Calculate slope_angle_max
!        In the x-direction
         do_it_x_u = .false.
         do_it_x_d = .false.
         number_pts = 0
         do iy=start_y+3,end_y-3
            do ix=start_x+3,end_x-3
               do iz=start_z,end_z
        	  slope_angle_x(ix,iy,iz) = abs(array(ix,iy,iz) - array(ix-1,iy,iz)) / resolution_x
        	  slope_angle_x(ix,iy,iz) = atan(slope_angle_x(ix,iy,iz)) * 180. / (4 * atan(1.))
        	  if (slope_angle_x(ix,iy,iz) > max_slope_angle(id-1) + epsilon_slope_angle) then 
        	     number_pts = number_pts + 1
        	     if (array(ix,iy,iz) > array(ix-1,iy,iz)) then
        		do_it_x_u(ix,iy,iz) = .false.
        		do_it_x_d(ix,iy,iz) = .true.
        	     else
        		do_it_x_u(ix,iy,iz) = .true.
        		do_it_x_d(ix,iy,iz) = .false.
        	     end if
        	  end if
               end do
            end do
         end do
         call mprintf(.true.,LOGFILE,'    %i point(s) exceed(s) the target angle (%f) in the x-direction',i1=number_pts,f1=max_slope_angle(id-1))
         if ((number_pts > 0) .and. (debug_level >= 100)) then
     	    do iy=start_y+3,end_y-3
     	       do ix=start_x+3,end_x-3
     	  	  do iz=start_z,end_z
     	  	     if (do_it_x_u(ix,iy,iz) .or. do_it_x_d(ix,iy,iz)) then
     	  		call mprintf(.true.,LOGFILE,'         at point: i=%i, j=%i',i1=ix+3,i2=iy+3)
     	  	     end if
     	  	  end do
     	       end do
     	    end do            
         end if
!        In the y-direction
         do_it_y_u = .false.
         do_it_y_d = .false.
         number_pts = 0
         do iy=start_y+3,end_y-3
            do ix=start_x+3,end_x-3
               do iz=start_z,end_z
        	  slope_angle_y(ix,iy,iz) = abs(array(ix,iy,iz) - array(ix,iy-1,iz)) / resolution_y
        	  slope_angle_y(ix,iy,iz) = atan(slope_angle_y(ix,iy,iz)) * 180. / (4 * atan(1.))
        	  if (slope_angle_y(ix,iy,iz) > max_slope_angle(id-1) + epsilon_slope_angle) then
        	     number_pts = number_pts + 1
        	     if (array(ix,iy,iz) > array(ix,iy-1,iz)) then
        		do_it_y_u(ix,iy,iz) = .false.
        		do_it_y_d(ix,iy,iz) = .true.
        	     else
        		do_it_y_u(ix,iy,iz) = .true.
        		do_it_y_d(ix,iy,iz) = .false.
        	     end if
        	  end if
               end do
            end do
         end do
         call mprintf(.true.,LOGFILE,'    %i point(s) exceed(s) the target angle (%f) in the y-direction',i1=number_pts,f1=max_slope_angle(id-1))
         if ((number_pts > 0) .and. (debug_level >= 100)) then
     	    do iy=start_y+3,end_y-3
     	       do ix=start_x+3,end_x-3
     	  	  do iz=start_z,end_z
     	  	     if (do_it_y_u(ix,iy,iz) .or. do_it_y_d(ix,iy,iz)) then
     	  		call mprintf(.true.,LOGFILE,'         at point: i=%i, j=%i',i1=ix+3,i2=iy+3)
     	  	     end if
     	  	  end do
     	       end do
     	    end do            
         end if
     	 number_pts = 0
     	 do iy=start_y+3,end_y-3
     	    do ix=start_x+3,end_x-3
     	       do iz=start_z,end_z
     		  if ((do_it_x_d(ix,iy,iz) .and. do_it_y_u(ix,iy,iz)) .or. &
     		      (do_it_x_u(ix,iy,iz) .and. do_it_y_d(ix,iy,iz))) then
     		     number_pts = number_pts + 1
     		  end if
     	       end do
     	    end do
     	 end do
     	 call mprintf(.true.,LOGFILE,'    %i saddle point(s) exceed(s) the target angle (%f) in both the x- and y-directions',i1=number_pts,f1=max_slope_angle(id-1))
         if ((number_pts > 0) .and. (debug_level >= 100)) then
     	    do iy=start_y+3,end_y-3
     	       do ix=start_x+3,end_x-3
     	  	  do iz=start_z,end_z
     	  	     if ((do_it_x_d(ix,iy,iz) .and. do_it_y_u(ix,iy,iz)) .or. &
     	  		 (do_it_x_u(ix,iy,iz) .and. do_it_y_d(ix,iy,iz))) then
     	  		call mprintf(.true.,LOGFILE,'         at point: i=%i, j=%i',i1=ix+3,i2=iy+3)
     	  	     end if
     	  	  end do
     	       end do
     	    end do            
         end if
         slope_angle_x_max = MAXVAL(slope_angle_x(start_x+3:end_x-3,start_y+3:end_y-3,:))
         slope_angle_y_max = MAXVAL(slope_angle_y(start_x+3:end_x-3,start_y+3:end_y-3,:))
         slope_angle_max = MAX(slope_angle_x_max,slope_angle_y_max)

         call mprintf(.true.,LOGFILE,'    Maximum angle after scaling: %f',f1=slope_angle_max)

      end if

      counter = 1

      do while (slope_angle_max > (max_slope_angle(id-1) + epsilon_slope_angle) .and. (counter <= max_counter))

         call mprintf(.true.,LOGFILE,'    Begin iteration no. %i',i1=counter)

!        Smooth saddle points, that is
!        smooth the topography locally if (do_it_x_d(ix,iy,iz) .and. do_it_y_u(ix,iy,iz)) .or.
!                                         (do_it_x_u(ix,iy,iz) .and. do_it_y_d(ix,iy,iz)) = .true.
!        using a 9-point filter with the coefficient of the central point equal to 0.25
         call mprintf(.true.,LOGFILE,'         Smooth saddle points using a 9-point filter')
         do iy=start_y+2,end_y-2
            do ix=start_x+2,end_x-2
               do iz=start_z,end_z
             	  if ((do_it_x_d(ix,iy,iz) .and. do_it_y_u(ix,iy,iz)) .or. &
                      (do_it_x_u(ix,iy,iz) .and. do_it_y_d(ix,iy,iz))) then
                     array(ix,iy,iz) = array(ix,iy,iz) / 4. &
                                     + (array(ix-1,iy,iz) + array(ix+1,iy,iz) + array(ix,iy-1,iz) + array(ix,iy+1,iz)) / 8. &
                                     + (array(ix-1,iy-1,iz) + array(ix+1,iy-1,iz) + array(ix-1,iy+1,iz) + array(ix+1,iy+1,iz)) / 16.
             	  end if
               end do
            end do
         end do

!        Reduce any angle above the specified threshold
!        In the x-direction
         call mprintf(.true.,LOGFILE,'         Reduce any angle above the threshold in the x-direction')
         do iy=start_y+2,end_y-2
            do ix=start_x+2,end_x-2
               do iz=start_z,end_z
             	  slope_angle_x(ix,iy,iz) = abs(array(ix,iy,iz) - array(ix-1,iy,iz)) / resolution_x
             	  slope_angle_x(ix,iy,iz) = atan(slope_angle_x(ix,iy,iz)) * 180. / (4 * atan(1.))
             	  if (slope_angle_x(ix,iy,iz) > max_slope_angle(id-1)) then
             	     do_it_x(ix,iy,iz) = .true.
             	  else
             	     do_it_x(ix,iy,iz) = .false.
             	  end if
               end do
            end do
         end do
         do iy=start_y+2,end_y-2
            do ix=start_x+2,end_x-2
               do iz=start_z,end_z
                  if (do_it_x(ix,iy,iz)) then
                     if (array(ix,iy,iz) > array(ix-1,iy,iz)) then
                  	scratch(ix,iy,iz) = array(ix-1,iy,iz) + resolution_x * tan(max_slope_angle(id-1) * (4 * atan(1.)) / 180.)
                     else
                  	scratch(ix,iy,iz) = array(ix-1,iy,iz) - resolution_x * tan(max_slope_angle(id-1) * (4 * atan(1.)) / 180.)    
                     end if
                  else
                     scratch(ix,iy,iz) = array(ix,iy,iz)
                  end if
               end do
            end do
         end do
!        In the y-direction
         call mprintf(.true.,LOGFILE,'         Reduce any angle above the threshold in the y-direction')
         do iy=start_y+3,end_y-3
            do ix=start_x+3,end_x-3
               do iz=start_z,end_z
                  slope_angle_y(ix,iy,iz) = abs(scratch(ix,iy,iz) - scratch(ix,iy-1,iz)) / resolution_y
                  slope_angle_y(ix,iy,iz) = atan(slope_angle_y(ix,iy,iz)) * 180. / (4 * atan(1.))
                  if (slope_angle_y(ix,iy,iz) > max_slope_angle(id-1)) then
                     do_it_y(ix,iy,iz) = .true.
                  else
                     do_it_y(ix,iy,iz) = .false.
                  end if
               end do
            end do
         end do
         do iy=start_y+3,end_y-3
            do ix=start_x+3,end_x-3
               do iz=start_z,end_z
                  if (do_it_y(ix,iy,iz)) then
                     if (scratch(ix,iy,iz) > scratch(ix,iy-1,iz)) then
                  	array(ix,iy,iz) = scratch(ix,iy-1,iz) + resolution_y * tan(max_slope_angle(id-1) * (4 * atan(1.)) / 180.)
                     else
                  	array(ix,iy,iz) = scratch(ix,iy-1,iz) - resolution_y * tan(max_slope_angle(id-1) * (4 * atan(1.)) / 180.)      
                     end if
                  else
                     array(ix,iy,iz) = scratch(ix,iy,iz)
                  end if
               end do
            end do
         end do

 	 call exchange_halo_r(array, &
 			      start_x, end_x, start_y, end_y, start_z, end_z, &
 			      start_dom_x, end_dom_x, start_dom_y, end_dom_y, start_z, end_z)

         call mprintf(.true.,LOGFILE,'    End iteration no. %i',i1=counter)

         counter = counter + 1

!        Calculate slope_angle_max
!        In the x-direction
         do_it_x_u = .false.
         do_it_x_d = .false.
         number_pts = 0
         do iy=start_y+3,end_y-3
            do ix=start_x+3,end_x-3
               do iz=start_z,end_z
        	  slope_angle_x(ix,iy,iz) = abs(array(ix,iy,iz) - array(ix-1,iy,iz)) / resolution_x
        	  slope_angle_x(ix,iy,iz) = atan(slope_angle_x(ix,iy,iz)) * 180. / (4 * atan(1.))
        	  if (slope_angle_x(ix,iy,iz) > max_slope_angle(id-1) + epsilon_slope_angle) then 
        	     number_pts = number_pts + 1
        	     if (array(ix,iy,iz) > array(ix-1,iy,iz)) then
        		do_it_x_u(ix,iy,iz) = .false.
        		do_it_x_d(ix,iy,iz) = .true.
        	     else
        		do_it_x_u(ix,iy,iz) = .true.
        		do_it_x_d(ix,iy,iz) = .false.
        	     end if
        	  end if
               end do
            end do
         end do
         call mprintf(.true.,LOGFILE,'    %i point(s) exceed(s) the target angle (%f) in the x-direction',i1=number_pts,f1=max_slope_angle(id-1))
         if ((number_pts > 0) .and. (debug_level >= 100)) then
     	    do iy=start_y+3,end_y-3
     	       do ix=start_x+3,end_x-3
     	  	  do iz=start_z,end_z
     	  	     if (do_it_x_u(ix,iy,iz) .or. do_it_x_d(ix,iy,iz)) then
     	  		call mprintf(.true.,LOGFILE,'         at point: i=%i, j=%i',i1=ix+3,i2=iy+3)
     	  	     end if
     	  	  end do
     	       end do
     	    end do            
         end if
!        In the y-direction
         do_it_y_u = .false.
         do_it_y_d = .false.
         number_pts = 0
         do iy=start_y+3,end_y-3
            do ix=start_x+3,end_x-3
               do iz=start_z,end_z
        	  slope_angle_y(ix,iy,iz) = abs(array(ix,iy,iz) - array(ix,iy-1,iz)) / resolution_y
        	  slope_angle_y(ix,iy,iz) = atan(slope_angle_y(ix,iy,iz)) * 180. / (4 * atan(1.))
        	  if (slope_angle_y(ix,iy,iz) > max_slope_angle(id-1) + epsilon_slope_angle) then
        	     number_pts = number_pts + 1
        	     if (array(ix,iy,iz) > array(ix,iy-1,iz)) then
        		do_it_y_u(ix,iy,iz) = .false.
        		do_it_y_d(ix,iy,iz) = .true.
        	     else
        		do_it_y_u(ix,iy,iz) = .true.
        		do_it_y_d(ix,iy,iz) = .false.
        	     end if
        	  end if
               end do
            end do
         end do
         call mprintf(.true.,LOGFILE,'    %i point(s) exceed(s) the target angle (%f) in the y-direction',i1=number_pts,f1=max_slope_angle(id-1))
         if ((number_pts > 0) .and. (debug_level >= 100)) then
     	    do iy=start_y+3,end_y-3
     	       do ix=start_x+3,end_x-3
     	  	  do iz=start_z,end_z
     	  	     if (do_it_y_u(ix,iy,iz) .or. do_it_y_d(ix,iy,iz)) then
     	  		call mprintf(.true.,LOGFILE,'         at point: i=%i, j=%i',i1=ix+3,i2=iy+3)
     	  	     end if
     	  	  end do
     	       end do
     	    end do            
         end if
     	 number_pts = 0
     	 do iy=start_y+3,end_y-3
     	    do ix=start_x+3,end_x-3
     	       do iz=start_z,end_z
     		  if ((do_it_x_d(ix,iy,iz) .and. do_it_y_u(ix,iy,iz)) .or. &
     		      (do_it_x_u(ix,iy,iz) .and. do_it_y_d(ix,iy,iz))) then
     		     number_pts = number_pts + 1
     		  end if
     	       end do
     	    end do
     	 end do
     	 call mprintf(.true.,LOGFILE,'    %i saddle point(s) exceed(s) the target angle (%f) in both the x- and y-directions',i1=number_pts,f1=max_slope_angle(id-1))
         if ((number_pts > 0) .and. (debug_level >= 100)) then
     	    do iy=start_y+3,end_y-3
     	       do ix=start_x+3,end_x-3
     	  	  do iz=start_z,end_z
     	  	     if ((do_it_x_d(ix,iy,iz) .and. do_it_y_u(ix,iy,iz)) .or. &
     	  		 (do_it_x_u(ix,iy,iz) .and. do_it_y_d(ix,iy,iz))) then
     	  		call mprintf(.true.,LOGFILE,'         at point: i=%i, j=%i',i1=ix+3,i2=iy+3)
     	  	     end if
     	  	  end do
     	       end do
     	    end do            
         end if
         slope_angle_x_max = MAXVAL(slope_angle_x(start_x+3:end_x-3,start_y+3:end_y-3,:))
         slope_angle_y_max = MAXVAL(slope_angle_y(start_x+3:end_x-3,start_y+3:end_y-3,:))
         slope_angle_max = MAX(slope_angle_x_max,slope_angle_y_max)

         call mprintf(.true.,LOGFILE,'    Maximum angle: %f',f1=slope_angle_max)

      end do ! do while loop

      if ((counter-1) == max_counter) then
         call mprintf(.true.,LOGFILE,'    Maximum number of iterations reached')
      end if

!     Smooth the topography locally for any angle within 10% of the specified threshold
!     using a 9-point filter with the coefficient of the central point equal to 0.25


      call mprintf(.true.,LOGFILE,'    Round any angle within 10 percents of the specified threshold using a 9-point filter')

      number_pts = 0

!     In the x-direction
      if (debug_level >= 100) then
         call mprintf(.true.,LOGFILE,'    In the x-direction')
      end if
      do iy=start_y+3,end_y-3
         do ix=start_x+3,end_x-3
            do iz=start_z,end_z
               slope_angle_x(ix,iy,iz) = abs(array(ix,iy,iz) - array(ix-1,iy,iz)) / resolution_x
               slope_angle_x(ix,iy,iz) = atan(slope_angle_x(ix,iy,iz)) * 180. / (4 * atan(1.))
               if (slope_angle_x(ix,iy,iz) > (max_slope_angle(id-1) * 0.9)) then
                  number_pts = number_pts + 1
                  if (debug_level >= 100) then
     	             call mprintf(.true.,LOGFILE,'         at point: i=%i, j=%i',i1=ix+3,i2=iy+3)
                  end if
                  array(ix,iy,iz) = array(ix,iy,iz) / 4. &
                  		  + (array(ix-1,iy,iz) + array(ix+1,iy,iz) + array(ix,iy-1,iz) + array(ix,iy+1,iz)) / 8. &
                  		  + (array(ix-1,iy-1,iz) + array(ix+1,iy-1,iz) + array(ix-1,iy+1,iz) + array(ix+1,iy+1,iz)) / 16.
               end if
            end do
         end do
      end do

!     In the y-direction
      if (debug_level >= 100) then
         call mprintf(.true.,LOGFILE,'    In the y-direction')
      end if
      do iy=start_y+3,end_y-3
         do ix=start_x+3,end_x-3
            do iz=start_z,end_z
               slope_angle_y(ix,iy,iz) = abs(scratch(ix,iy,iz) - scratch(ix,iy-1,iz)) / resolution_y
               slope_angle_y(ix,iy,iz) = atan(slope_angle_y(ix,iy,iz)) * 180. / (4 * atan(1.))
               if (slope_angle_y(ix,iy,iz) > (max_slope_angle(id-1) * 0.9)) then
                  number_pts = number_pts + 1
                  if (debug_level >= 100) then
     	             call mprintf(.true.,LOGFILE,'         at point: i=%i, j=%i',i1=ix+3,i2=iy+3)
                  end if
                  array(ix,iy,iz) = array(ix,iy,iz) / 4. &
                  		  + (array(ix-1,iy,iz) + array(ix+1,iy,iz) + array(ix,iy-1,iz) + array(ix,iy+1,iz)) / 8. &
                  		  + (array(ix-1,iy-1,iz) + array(ix+1,iy-1,iz) + array(ix-1,iy+1,iz) + array(ix+1,iy+1,iz)) / 16.
               end if
            end do
         end do
      end do
      call mprintf(.true.,LOGFILE,'    Total number of point(s): %i',i1=number_pts)

      call exchange_halo_r(array, &
 	     		   start_x, end_x, start_y, end_y, start_z, end_z, &
 	     		   start_dom_x, end_dom_x, start_dom_y, end_dom_y, start_z, end_z)


!     Remove artificially negative values
      do iy=start_y,end_y
 	 do ix=start_x,end_x
 	    do iz=start_z,end_z
 	       if (array(ix,iy,iz) < 0. .and. orig_array(ix,iy,iz) >= 0.) then
 		  array(ix,iy,iz) = orig_array(ix,iy,iz)
 	       end if
 	    end do
 	 end do
      end do

!     Apply npass smth_opt passes
      if (smth_opt == ONETWOONE) then
 	 call one_two_one(array,		  &
 			  start_dom_x, end_dom_x, &
 			  start_dom_y, end_dom_y, &
 			  start_x,     end_x,	  &
 			  start_y,     end_y,	  &
 			  start_z,     end_z,	  &
 			  npass, msgval)
      else if (smth_opt == SMTHDESMTH) then
 	 call smth_desmth(array,		  &
 			  start_dom_x, end_dom_x, &
 			  start_dom_y, end_dom_y, &
 			  start_x,     end_x,	  &
 			  start_y,     end_y,	  &
 			  start_z,     end_z,	  &
 			  npass, msgval)
      else if (smth_opt == SMTHDESMTH_SPECIAL) then
 	 call smth_desmth_special(array,	  &
 			  start_dom_x, end_dom_x, &
 			  start_dom_y, end_dom_y, &
 			  start_x,     end_x,	  &
 			  start_y,     end_y,	  &
 			  start_z,     end_z,	  &
 			  npass, msgval)
      end if

      deallocate(scratch)
      deallocate(orig_array)
!     Deallocate new variables
      deallocate(slope_angle_x)
      deallocate(slope_angle_y)
      deallocate(do_it_x)
      deallocate(do_it_x_u)
      deallocate(do_it_x_d)
      deallocate(do_it_y)
      deallocate(do_it_y_u)
      deallocate(do_it_y_d)

   end subroutine smth_cc_ardu
!E: Added by CC

   !
   ! Smoothing routines for E-grid, contributed by Matthew Pyle
   !

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! Name: one_two_one_egrid
   !
   ! Purpose: Apply the 1-2-1 smoother from the MM5 program TERRAIN 
   !   (found in smth121.F) to array.
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine one_two_one_egrid(array, start_dom_x, end_dom_x, start_dom_y, end_dom_y, &
                                start_x, end_x, start_y, end_y, start_z, end_z, npass, msgval, hflag)

      implicit none

      ! Arguments
      integer, intent(in) :: start_dom_x, start_dom_y, start_x, start_y, start_z
      integer, intent(in) :: end_dom_x, end_dom_y, end_x, end_y, end_z
      integer, intent(in) :: npass
      real, intent(in) :: msgval, hflag
      real, dimension(start_x:end_x, start_y:end_y, start_z:end_z), intent(inout) :: array

      ! Local variables
      integer :: ix, iy, iz, ipass
      real, pointer, dimension(:,:,:) :: scratch
      integer, dimension(start_y:end_y) :: ihe, ihw, istart

      allocate(scratch(start_x:end_x, start_y:end_y, start_z:end_z))

      do iy=start_y,end_y
         if (hflag == 1.0) then
            ihe(iy) = abs(mod(iy+1,2))
            ihw(iy) = ihe(iy)-1
         else
            ! assign ive,ivw equivs to ihe,ihw
            ihe(iy) = abs(mod(iy,2))
            ihw(iy) = ihe(iy)-1
         end if
      end do

      do iy=start_y,end_y
         if (hflag == 1.0) then
            if (mod(iy,2) == 0) then
               istart(iy) = start_x
            else
               istart(iy) = start_x+1
            end if
         else ! v points
            if (abs(mod(iy,2)) == 1) then
               istart(iy) = start_x
            else
               istart(iy) = start_x+1
            end if
         end if
      end do

      do ipass=1,npass

         do iy=start_y,end_y
            do ix=start_x,end_x
               scratch(ix,iy,1) = array(ix,iy,1) ! for points used in 2nd computation but not defined in 1st computation
            end do
         end do

         ! SW-NE direction
         do iy=start_y+1,end_y-1
            do ix=istart(iy),end_x-1
               do iz=start_z,end_z
                  if ( (msgval == 1.0 .and. array(ix,iy,iz) /= 0.) .or. msgval /= 1.0) then
                     scratch(ix,iy,iz) = 0.50*array(ix,iy,iz)+ &
                                      0.25*(array(ix+ihw(iy),iy-1,iz)+array(ix+ihe(iy),iy+1,iz))
                  end if
               end do
            end do
         end do

         ! NW-SE direction
         do iy=start_y+1,end_y-1
            do ix=istart(iy),end_x-1
               do iz=start_z,end_z
                  if ( (msgval == 1.0 .and. array(ix,iy,iz) /= 0.) .or. msgval /= 1.0) then
                     array(ix,iy,iz) = 0.50*scratch(ix,iy,iz)+ &
                                    0.25*(scratch(ix+ihe(iy),iy-1,iz)+scratch(ix+ihw(iy),iy+1,iz))
                  end if
               end do
            end do
         end do

         call exchange_halo_r(array, &
                              start_x, end_x, start_y, end_y, start_z, end_z, &
                              start_dom_x, end_dom_x, start_dom_y, end_dom_y, start_z, end_z)

      end do

      deallocate(scratch)

   end subroutine one_two_one_egrid


   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! Name: smth_desmth_egrid_old
   !
   ! Purpose: Apply the smoother-desmoother for E grid
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine smth_desmth_egrid_old(array, start_dom_x, end_dom_x, start_dom_y, end_dom_y, &
                                    start_x, end_x, start_y, end_y, start_z, end_z, npass, msgval, hflag)

      implicit none

      ! Arguments
      integer, intent(in) :: start_dom_x, start_dom_y, start_x, start_y, start_z
      integer, intent(in) :: end_dom_x, end_dom_y, end_x, end_y, end_z
      integer, intent(in) :: npass
      real, intent(in) :: msgval, hflag
      real, dimension(start_x:end_x, start_y:end_y, start_z:end_z), &
               intent(inout) :: array

      ! Local variables
      integer :: ix, iy, iz, ipass
      real, pointer, dimension(:,:,:) :: scratch
      integer, dimension(start_y:end_y) :: ihe, ihw, istart
      real, parameter:: cenwgt = 1.52
      real, parameter:: endwgt = 0.13

      allocate(scratch(start_x:end_x, start_y:end_y, start_z:end_z))

      do iy=start_y,end_y
         if (hflag == 1.0) then
            ihe(iy) = abs(mod(iy+1,2))
            ihw(iy) = ihe(iy)-1
         else
            ! assign ive,ivw equivs to ihe,ihw
            ihe(iy) = abs(mod(iy,2))
            ihw(iy) = ihe(iy)-1
         end if
      end do

      do iy=start_y,end_y
         if (hflag == 1.0) then
            if (mod(iy,2) == 0) then
               istart(iy) = start_x
            else
               istart(iy) = start_x+1
            endif
         else ! v points
            if (abs(mod(iy,2)) == 1) then
               istart(iy) = start_x
            else
               istart(iy) = start_x+1
            endif
         endif
      end do

      do ipass=1,npass

         !
         ! Smoothing pass
         !

         do iy=start_y,end_y
            do ix=start_x,end_x
               scratch(ix,iy,1) = array(ix,iy,1) 
            end do
         end do

         do iy=start_y+1,end_y-1
            do ix=istart(iy),end_x-1
               do iz=start_z,end_z
                  if ( (msgval == 1.0 .and. array(ix,iy,iz) /= 0.) .or. msgval /= 1.0) then
                     scratch(ix,iy,iz) = 0.50*array(ix,iy,iz)+ &
                                      0.125*(array(ix+ihw(iy),iy-1,iz)+array(ix+ihe(iy),iy+1,iz)+ &
                                             array(ix+ihw(iy),iy+1,iz)+array(ix+ihe(iy),iy-1,iz))
                  end if
               end do
            end do
         end do


         !
         ! Desmoothing pass
         !

         do iy=start_y+2,end_y-2
            do ix=istart(iy),end_x-1
               do iz=start_z,end_z
                  if ( (msgval == 1.0 .and. scratch(ix,iy,iz) /= 0.) .or. msgval /= 1.0) then
                     array(ix,iy,iz) = cenwgt*scratch(ix,iy,iz) - &
                                      endwgt*(scratch(ix+ihw(iy),iy-1,iz)+scratch(ix+ihe(iy),iy+1,iz) + &
                                              scratch(ix+ihw(iy),iy+1,iz)+scratch(ix+ihe(iy),iy-1,iz))
                  end if
               end do
            end do
         end do

      end do

      deallocate(scratch)

   end subroutine smth_desmth_egrid_old


   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! Name: smth_desmth_egrid
   !
   ! Purpose: Apply the smoother-desmoother for E grid 
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine smth_desmth_egrid(array, start_dom_x, end_dom_x, start_dom_y, end_dom_y, &
                                start_x, end_x, start_y, end_y, start_z, end_z, npass, msgval, hflag)

      implicit none

      ! Arguments
      integer, intent(in) :: start_dom_x, start_dom_y, start_x, start_y, start_z
      integer, intent(in) :: end_dom_x, end_dom_y, end_x, end_y, end_z
      integer, intent(in) :: npass
      real, intent(in) :: msgval, hflag
      real, dimension(start_x:end_x, start_y:end_y, start_z:end_z), &
               intent(inout) :: array

      ! Local variables
      integer :: ix, iy, iz, ipass
      real, pointer, dimension(:,:,:) :: scratch
      integer, dimension(start_y:end_y) :: ihe, ihw, istart
      real, parameter :: cenwgt = 1.52
      real, parameter :: endwgt = 0.26

      allocate(scratch(start_x:end_x, start_y:end_y, start_z:end_z))

      do iy=start_y,end_y

         if (hflag .eq. 1.0) then
            ihe(iy)=abs(mod(iy+1,2))
            ihw(iy)=ihe(iy)-1

         ! assign ive,ivw equivs to ihe,ihw
         else
            ihe(iy)=abs(mod(iy,2))
            ihw(iy)=ihe(iy)-1

         end if

      end do

      do iy=start_y,end_y

         if (hflag .eq. 1.0) then
            if (mod(iy,2) .eq. 0) then
               istart(iy)=start_x
            else
               istart(iy)=start_x+1
            endif

         else ! v points
            if (abs(mod(iy,2)) .eq. 1) then
               istart(iy)=start_x
            else
               istart(iy)=start_x+1
            end if

         end if

      end do


      do ipass=1,npass

         !
         ! Smoothing pass
         !

         do iy=start_y,end_y
         do ix=start_x,end_x
            scratch(ix,iy,1)=array(ix,iy,1) ! for points used in 2nd computation but 
                                            !    not defined in 1st
         end do
         end do

         ! SW-NE direction
         do iy=start_y+1,end_y-1
            do ix=istart(iy),end_x-1
               do iz=start_z,end_z
                  if ( (msgval .eq. 1.0 .and. array(ix,iy,iz) .ne. 0.) .or. msgval .ne. 1.0) then
                     scratch(ix,iy,iz) = 0.50*array(ix,iy,iz)+ &
                     0.25*(array(ix+ihw(iy),iy-1,iz)+array(ix+ihe(iy),iy+1,iz))
                  end if
               end do
            end do
         end do

         ! NW-SE direction
         do iy=start_y+1,end_y-1
            do ix=istart(iy),end_x-1
               do iz=start_z,end_z
                  if ( (msgval .eq. 1.0 .and. array(ix,iy,iz) .ne. 0.) .or. msgval .ne. 1.0) then
                     array(ix,iy,iz) = 0.50*scratch(ix,iy,iz)+ &
                     0.25*(scratch(ix+ihe(iy),iy-1,iz)+scratch(ix+ihw(iy),iy+1,iz))
                  end if
               end do
            end do
         end do

         call exchange_halo_r(array, &
                              start_x, end_x, start_y, end_y, start_z, end_z, &
                              start_dom_x, end_dom_x, start_dom_y, end_dom_y, start_z, end_z)



         !
         ! Desmoothing pass
         !

         ! SW-NE direction
         do iy=start_y+2,end_y-2
            do ix=istart(iy),end_x-1
               do iz=start_z,end_z
                  if ( (msgval .eq. 1.0 .and. array(ix,iy,iz) .ne. 0.) .or. msgval .ne. 1.0) then
                     scratch(ix,iy,iz) = cenwgt*array(ix,iy,iz) - &
                       endwgt*(array(ix+ihw(iy),iy-1,iz)+array(ix+ihe(iy),iy+1,iz))
                  end if
               end do
            end do
         end do

         ! NW-SE direction
         do iy=start_y+2,end_y-2
            do ix=istart(iy),end_x-1
               do iz=start_z,end_z
                  if ( (msgval .eq. 1.0 .and. array(ix,iy,iz) .ne. 0.) .or. msgval .ne. 1.0) then
                     array(ix,iy,iz) = cenwgt*scratch(ix,iy,iz) - &
                       endwgt*(scratch(ix+ihe(iy),iy-1,iz)+scratch(ix+ihw(iy),iy+1,iz))
                  end if
               end do
            end do
         end do

         call exchange_halo_r(array, &
                              start_x, end_x, start_y, end_y, start_z, end_z, &
                              start_dom_x, end_dom_x, start_dom_y, end_dom_y, start_z, end_z)

      end do

      deallocate(scratch)

   end subroutine smth_desmth_egrid

end module smooth_module

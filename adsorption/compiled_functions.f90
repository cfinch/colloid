module compiled_functions
! Build this module with the command:
!     f2py -c -m compiled compiled_functions.f90
! With Intel compiler:
!     f2py --fcompiler=intelem -c -m compiled compiled_functions.f90 

      contains

      function compute_avf(grid_len, x_grid, y_grid, h_len, h_grid, &
        num_adsorbed, adsorbed_x, adsorbed_y, r, num_replicates, random_radius)
            implicit none

            integer, intent(in) :: grid_len, num_adsorbed, h_len, num_replicates
            real, intent(in), dimension(grid_len) :: x_grid, y_grid
            real, intent(in), dimension(h_len) :: h_grid
            real, intent(in) :: r, random_radius
            real, intent(in), dimension(num_adsorbed) :: adsorbed_x, adsorbed_y
            !f2py depend(grid_len) x_grid, y_grid            
            !f2py depend(h_len) h_grid, compute_avf
            !f2py depend(num_adsorbed) adsorbed_x, adsorbed_y

            real, dimension(h_len) :: compute_avf

            ! Local variables
            integer :: hstep, replicate, i,j,p      ! Counters
            integer :: num_available
            real :: h,r_squared,rand
            real, dimension(grid_len, grid_len) :: x,y
            real, dimension(h_len) :: avf_data

            r_squared = 4.0*r**2
            avf_data(:) = 0.0
            
            do replicate = 1,num_replicates
                ! Generate new randomized positions for test particles
                do i=1,grid_len
                    do j=1,grid_len
                        call random_number(rand)
                        x(i,j) = x_grid(i) + (rand-0.5)*random_radius
                        if (x(i,j) > x_grid(grid_len)) then
                            x(i,j) = x_grid(grid_len)
                        else if (x(i,j) < x_grid(1)) then
                            x(i,j) = x_grid(1)
                        endif
                        call random_number(rand)
                        y(i,j) = y_grid(j) + (rand-0.5)*random_radius
                        if (y(i,j) > y_grid(grid_len)) then
                            y(i,j) = y_grid(grid_len)
                        else if (y(i,j) < y_grid(1)) then
                            y(i,j) = y_grid(1)
                        endif
                    enddo
                enddo

                ! Use these test particles to compute AVF at all h locations
                do hstep = 1,h_len
                    h = h_grid(hstep)
                    num_available = grid_len**2
                    do i=1,grid_len
                        do j=1,grid_len
                            do p=1,num_adsorbed
                                if ((x(i,j)-adsorbed_x(p))**2 + (y(i,j)-adsorbed_y(p))**2 + h**2 < r_squared) then
                                    num_available = num_available - 1
                                    exit
                                endif
                            enddo   ! p
                        enddo   ! j
                    enddo   ! i
                    avf_data(hstep) = avf_data(hstep) + real(num_available)/real(grid_len**2)
                enddo   ! hstep
            enddo   ! replicate
           
            forall (hstep=1:h_len) compute_avf(hstep) = avf_data(hstep)/real(num_replicates)

      end function compute_avf

end module compiled_functions


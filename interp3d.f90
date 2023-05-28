! interp3d.f90 
! Interpolate EM density grids

! Copyright (C) 2021 Jiahua He, Tao Li, Sheng-You Huang and Huazhong University of Science and Technology
! This program is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.

! This program is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.

! You should have received a copy of the GNU General Public License
! along with this program.  If not, see <http://www.gnu.org/licenses/>.

module interp3d
real(kind=4), allocatable :: mapout(:, :, :)
integer pextx, pexty, pextz
contains
    subroutine linear(mapin, zpix, ypix, xpix, apix, nz, ny, nx)
    implicit none
    real zpix, ypix, xpix, apix
    integer nx, ny, nz
    real mapin(nz, ny, nx)
    integer indx, indy, indz
    real xpos, ypos, zpos, gx, gy, gz, a, b, c
    integer x0, y0, z0, x1, y1, z1

    pextx = floor(xpix * (nx - 1) / apix) + 1
    pexty = floor(ypix * (ny - 1) / apix) + 1
    pextz = floor(zpix * (nz - 1) / apix) + 1
    allocate(mapout(pextz, pexty, pextx))
    do indz = 1, pextz
        do indy = 1, pexty
            do indx = 1, pextx
                xpos = (indx - 1) * apix
                ypos = (indy - 1) * apix
                zpos = (indz - 1) * apix
                gx = xpos / xpix + 1
                gy = ypos / ypix + 1
                gz = zpos / zpix + 1
                x0 = floor(gx)
                y0 = floor(gy)
                z0 = floor(gz)
                x1 = x0 + 1
                y1 = y0 + 1
                z1 = z0 + 1
                if(x0 >= 1 .and. x1 <= nx &
             .and. y0 >= 1 .and. y1 <= ny &
             .and. z0 >= 1 .and. z1 <= nz) then
                    a = gx - x0
                    b = gy - y0
                    c = gz - z0
           mapout(indz, indy, indx)=a    *b    *c    *mapin(z1,y1,x1) &
                                   +(1-a)*b    *c    *mapin(z1,y1,x0) &
                                   +a    *(1-b)*c    *mapin(z1,y0,x1) &
                                   +a    *b    *(1-c)*mapin(z0,y1,x1) &
                                   +a    *(1-b)*(1-c)*mapin(z0,y0,x1) &
                                   +(1-a)*b    *(1-c)*mapin(z0,y1,x0) &
                                   +(1-a)*(1-b)*c    *mapin(z1,y0,x0) &
                                   +(1-a)*(1-b)*(1-c)*mapin(z0,y0,x0)
                end if
            end do 
        end do 
    end do

    return
    end subroutine

    subroutine del_mapout
    if (allocated(mapout)) then
        deallocate(mapout)
    end if
    return 
    end subroutine

end module  

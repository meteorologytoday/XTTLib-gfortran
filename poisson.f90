
module poisson

  use xttdft

contains

  ! used with "poisson_inverse_fft"
  ! to avoid the recalculation of tmp
  subroutine make_poisson_mtx(P2mtx, dmn, Lx, Ly)
    integer :: dmn(2)
    complex :: P2mtx(dmn(1),dmn(2))
    real    :: Lx, Ly

    real    :: tmpx(dmn(1)), tmpy(dmn(2)), tmp

    ! for x
    if(mod(dmn(1), 2) .eq. 0) then
       ! even
       tmpx(1) = 0
       tmpx(dmn(1)/2+1) = 0
       
       do i=1, dmn(1)/2-1
          tmp = (i*PII/Lx)**2
          tmpx(1+i)    = tmp
          tmpx(dmn(1)+1-i) = tmp
       end do
       
    else
       ! odd
       tmpx(1) = 0
       do i=1, (dmn(1)-1)/2
          tmp = (i*PII/Lx)**2
          tmpx(1+i)    = tmp
          tmpx(dmn(1)+1-i) = tmp
       end do
       
    end if

    ! for y
    if(mod(dmn(2), 2) .eq. 0) then
       tmpy(1) = 0
       tmpy(dmn(2)/2+1) = 0
       
       do i=1, dmn(2)/2-1
          tmp = (i*PII/Ly)**2
          tmpy(1+i)    = tmp
          tmpy(dmn(2)+1-i) = tmp
       end do
       
    else
       tmpy(1) = 0
       do i=1, (dmn(2)-1)/2
          tmp = (i*PII/Ly)**2
          tmpy(1+i)    = tmp
          tmpy(dmn(2)+1-i) = tmp
       end do
    end if

 
    do i=1,dmn(1)
       do j=1,dmn(2)
          P2mtx(i,j) = - (tmpx(i) + tmpy(j))
       end do
    end do


  end subroutine make_poisson_mtx

  ! p   : poisson
  subroutine solve_poisson_dft(X, p, dmn, P2mtx, W2mtx, iW2mtx)
    integer  :: dmn(2)
    real     :: X(dmn(1),dmn(2)), p(dmn(1),dmn(2)), Lx, Ly


    complex  :: XX(dmn(1), dmn(2)), P2mtx(dmn(1),dmn(2)), W2mtx(dmn(1),dmn(2)), iW2mtx(dmn(1),dmn(2))
    integer  :: i,j
    
    call dft2(p,XX, dmn, W2mtx)

    ! do inverse


    do i=1,dmn(1)
       do j=1,dmn(2)
          if(P2mtx(i,j) .ne. 0) then
             XX(i,j) = XX(i,j) / P2mtx(i,j)
          else 
             XX(i,j) = 0
          end if
       end do
    end do

    ! idft2
    call idft2(X, XX, dmn, iW2mtx)
 
  end Subroutine solve_poisson_dft
  
  subroutine solve_poisson_jacobian_cyclic(f, g, dmn, delta, maxerr)
    integer :: dmn(2)
    real :: f(dmn(1), dmn(2)), f_next(dmn(1), dmn(2)), g(dmn(1), dmn(2)),  maxerr
    real :: g_num(dmn(1), dmn(2)), delta2
    
    delta2 = delta * delta
    
    call laplacian_diff_cyclic(g_num, f, delta)
    g_num = g_num - g
    
    f_next = f ! guess field
    
    do while (sum(g_num*g_num) .gt. maxerr)
       
       f_next = 0.25 * (cshift(f, shift = 1, dim = 1) + cshift(f, shift = -1, dim = 1) + &
                        cshift(f, shift = 1, dim = 2) + cshift(f, shift = -1, dim = 2) - delta2*g)
       
       f = f_next
       
       call laplacian_diff_cyclic(g_num, f, delta)
       g_num = g_num - g
       
    end do
    
    print *, "sum err square = ", sum(g_num*g_num)
  end subroutine solve_poisson_jacobian_cyclic
  
  
  
  subroutine laplacian_diff_cyclic(g, f, delta) 
    real :: f(:,:), g(:,:), delta
    
    g = (cshift(f, shift = 1, dim = 1) + cshift(f, shift = -1, dim = 1) + &
         cshift(f, shift = 1, dim = 2) + cshift(f, shift = -1, dim = 2) - 4*f) / (delta*delta)
  
  end subroutine laplacian_diff_cyclic

  subroutine solve_poisson_jacobian_noncyclic_well(f, g, dmn, delta, maxerr)
    integer :: dmn(2)
    real :: f(dmn(1), dmn(2)), f_next(dmn(1), dmn(2)), g(dmn(1), dmn(2)), filter(dmn(1), dmn(2)), maxerr
    real :: g_num(dmn(1), dmn(2)), delta2
    
    delta2 = delta * delta
    
    call laplacian_diff_cyclic(g_num, f, delta)
    g_num = g_num - g
    
    f_next = f ! guess field

    filter = 1
    filter(:,1) = 0
    filter(1,:) = 0
    filter(dmn(1), :) = 0
    filter(:, dmn(1)) = 0
    
    do while (sum(g_num*g_num) .gt. maxerr)
       
       f = f * filter

       f_next = 0.25 * (cshift(f, shift = 1, dim = 1) + cshift(f, shift = -1, dim = 1) + &
                        cshift(f, shift = 1, dim = 2) + cshift(f, shift = -1, dim = 2) - delta2*g)

       f = f_next
       
       call laplacian_diff_cyclic(g_num, f, delta)
       g_num = g_num - g
       
    end do
    
    print *, "sum err square = ", sum(g_num*g_num)
  end subroutine solve_poisson_jacobian_noncyclic_well

!  subroutine laplacian_diff_noncyclic(g, f, delta, dmn)
!    integer :: dmn(2)
!    real :: f(dmn(1), dmn(2)), g(dmn(1), dmn(2)), delta
    
!    g = (cshift(f, shift = 1, dim = 1) + cshift(f, shift = -1, dim = 1) + &
!         cshift(f, shift = 1, dim = 2) + cshift(f, shift = -1, dim = 2) - 4*f) / (delta*delta)
  
  !  g(1,:) = (g(3,:)+g(1,:) - 2*g(2,:))/(delta*delta)
 !   g(dmn(1),:) = (g(dmn(1))+)

!  end subroutine laplacian_diff_cyclic

 



!  ! FOR EQUAL SIZE(POWER OF 2) MATRIX
!  subroutine solve_poisson_jacobian_multigrid_cyclic(f, g, dmn, delta, maxerr)
!    integer :: dmn, subsize
!    real    :: f(dmn, dmn), g(dmn, dmn), f_fake(dmn, dmn), g_fake(dmn, dmn), delta, maxerr
!    real    :: avg

!  end subroutine solve_poisson_jacobian_multigrid_cyclic

!  subroutine solve_poisson_jacobian_multigrid_cyclic_sub(f, g, dmn, delta, maxerr)
    
!    subsize = dmn

!    if (subsize .ne. 1) then
!       subsize = subsize / 2

       ! cal mean
!       do i=1, dmn/subsize
!          do j=1, dmn/subsize
!             g_fake(i,j) = sum( g((i-1)*subsize+1:i*subsize, (j-1)*subsize+1:j*subsize)) / (subsize**2)
!          end do
!       end do

!       call solve_poisson_jacobian_cyclic(f_fake, g, dmn, delta, maxerr)

!    end do


!  end subroutine solve_poisson_jacobian_multigrid_cyclic_sub


end module poisson

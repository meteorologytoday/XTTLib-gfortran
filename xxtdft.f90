module xttdft

  real,    parameter :: PI = 3.141592653589793, PII = 6.283185307179586
  integer, parameter :: FOURIER_FORWARD = -1, FOURIER_BACKWARD = 1

contains

  subroutine dft2(X, XX, dmn, W2mtx)
    
    complex, parameter :: XJ = (0, FOURIER_FORWARD)

    integer            :: dmn(2)
    real               :: X(dmn(1), dmn(2))
    complex            :: XX(dmn(1), dmn(2)), W2mtx(dmn(1), dmn(2))
    
    integer            :: i, j, kx, ky, a, b
    integer            :: rows, cols
    
    
    rows = int((dmn(1)-1) / 2)
    cols = int((dmn(2)-1) / 2)
    
      
    do a=1, dmn(1)                          ! wave number space
       do b=1, dmn(2)-cols                  ! wave number space
          ! add all elements
          XX(a,b) = 0
          kx = 1
          ky = 1
          
          do i=1,dmn(1)
             do j=1,dmn(2)
                
                XX(a,b) = XX(a,b) + W2mtx(kx,ky) * X(i,j)
              
                ky = ky + (b-1) 
                
                if(ky > dmn(2)) then
                   ky = ky - dmn(2)
                end if
             end do
             kx = kx + (a-1)
             if(kx > dmn(1)) then
                kx = kx - dmn(1)
             end if
          end do
       end do
    end do
    
    
    ! mirror conj
    do b=2, cols+1
       XX(1,dmn(2)+2-b) = conjg(XX(1,b))
    end do
    
    do a=2, dmn(1)
       do b=2, cols+1
          ! 2 times mirror
          XX(dmn(1)+2-a,dmn(2)+2-b) = conjg(XX(a,b))
       end do
    end do
    
  end subroutine dft2
  
  
  subroutine idft2(X, XX, dmn, W2mtx)
    
    complex, parameter :: XJ = (0, FOURIER_BACKWARD)
    
    integer            :: dmn(2)
    real               :: X(dmn(1), dmn(2))
    complex            :: XX(dmn(1), dmn(2)), W2mtx(dmn(1), dmn(2))
    
    integer            :: i, j, kx, ky, a, b, N
    
    N = dmn(1) * dmn(2)
    rows = int((dmn(1)-1) / 2)
    cols = int((dmn(2)-1) / 2)
    
    
    do a=1, dmn(1)                          ! wave number space
       do b=1, dmn(2)                       ! wave number space
          ! add all elements
          X(a,b) = 0
          kx = 1
          ky = 1
          
          do i=1,dmn(1)
             do j=1,dmn(2)
                
                X(a,b) = X(a,b) + real(W2mtx(kx,ky) * XX(i,j))
                
                ky = ky + (b-1) 
                
                if(ky > dmn(2)) then
                   ky = ky - dmn(2)
                end if
             end do
             kx = kx + (a-1)
             if(kx > dmn(1)) then
                kx = kx - dmn(1)
             end if
          end do
          X(a,b) = X(a,b) / N
       end do
    end do
    
  end subroutine idft2
  
  subroutine make_Wmtx(exptable, dmn, dir)
    complex            :: exptable(dmn)
    integer            :: dmn, dir
    
    complex            :: XJ
    integer            :: i
    
    if (dir .eq. FOURIER_FORWARD) then
       XJ = (0,  FOURIER_FORWARD)
    else if (dir .eq. FOURIER_BACKWARD) then
       XJ = (0, FOURIER_BACKWARD)
    end if
    
    do i=1, dmn
       exptable(i) = exp(XJ*PII*(i-1)/dmn)
    end do
    
  end subroutine make_Wmtx
  
  
  subroutine make_W2mtx(exptable, dmn, dir)
    
    integer            :: dmn(2), dir
    complex            :: exptable(dmn(1), dmn(2))
    
    complex            :: XJ
    integer            :: i, j
    
    
    if (dir .eq. FOURIER_FORWARD) then
       XJ = (0,  FOURIER_FORWARD)
    else if (dir .eq. FOURIER_BACKWARD) then
       XJ = (0, FOURIER_BACKWARD)
    end if
    
    do i=1, dmn(1)
       do j=1, dmn(2)
          exptable(i,j) = exp(XJ*PII*(i-1)/dmn(1)) * exp(XJ*PII*(j-1)/dmn(2))
       end do
    end do
    
  end subroutine make_W2mtx
  
end module xttdft


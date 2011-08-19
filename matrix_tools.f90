module matrix_tools

  interface meshGrid
     module procedure meshGrid_2D
  end interface meshGrid

  interface clrMtx
     module procedure clrMtx_2D
  end interface clrMtx

  interface addGaussians
     module procedure addGaussians_2D
  end interface addGaussians

  interface addGaussian
     module procedure addGaussian_2D
  end interface addGaussian

contains




  function openLinspace(Xstart, Xend, N)
    real     :: openLinspace(N)
    real     :: Xstart, Xend, L
    integer  :: N,i

    L = Xend - Xstart
    N = abs(N)

    do i=1,N
       openLinspace(i) = Xstart + (L/N) * (i-1)
    end do
 
  end function openLinspace
  

  function closedLinspace(Xstart, Xend, N)
    real     :: closedLinspace(N)
    real     :: Xstart, Xend, L
    integer  :: N,i

    L = Xend - Xstart
    N = abs(N)

    if(N .eq. 1) then
       print *, 'Error, N must be grater than 1 !'
    else
       do i=1,N
          closedLinspace(i) = Xstart + (L/(N-1)) * (i-1)
       end do
    end if
  end function closedLinspace


  subroutine meshGrid_2D(X1, X2, X1pts, X2pts, Dim)
    real     :: X1(:,:), X2(:,:), X1pts(:), X2pts(:)
    integer  :: Dim(2), i, j

    do i=1,Dim(1)
       do j=1,Dim(2)
          X1(i,j) = X1pts(i)
          X2(i,j) = X2pts(j)
       end do
    end do
  end subroutine meshGrid_2D

  subroutine clrMtx_2D(data)
    real :: data(:,:)
    data = 0
  end subroutine clrMtx_2D
  
  ! X/Y        : the X/Ypts of the map
  ! Dim        : the rows and columns number of the map
  ! Gn         : the number of peaks
  ! CentX/Ys   : the centers of the gaussion
  ! Sigs       : the sigmas of the gaussion
  ! Maxs       : the maximas of the gaussion
  subroutine addGaussians_2D(data, X, Y, Dim, Gn, Cents, Sigs, Maxs)
    real          :: X(:,:), Y(:,:), data(:,:), Cents(:,:), Sigs(:), Maxs(:)
    integer       :: Dim(2) , Gn
    
    integer       :: n
    
    do n=1,Gn   
       call addGaussian_2D(data, X, Y, Dim, Cents(n,:), Sigs(n), Maxs(n))
    end do
    
  end subroutine addGaussians_2D
  
  subroutine addGaussian_2D(data, X1, X2, Dim, Cent, Sig, Max)
    real          :: X1(:,:), X2(:,:), data(:,:), Cent(2), Sig, Max
    integer       :: Dim(2)
    
    real          :: radius2
    integer       :: i,j
    
    do i=1, Dim(1)
       do j=1, Dim(2)
          radius2 = (X1(i,j) - Cent(1))**2 + (X2(i,j) - Cent(2))**2
          radius2 = - radius2 / (2 * Sig**2)
          data(i,j) = data(i,j) + Max * exp(radius2)
          print *, "(i,j) = (", i, ",",j, ") data = " , data(i,j), ", radius2 = " , radius2
       end do
    end do
    
  end subroutine addGaussian_2D
  
  
  subroutine mtxReadBinary_2D(data, unit, file, iostat)
    implicit none

    character(*) :: file
    real         :: data(:,:)
    integer      :: unit
  
    integer      :: iostat
    
    open(unit, file = trim(file), form = "unformatted", access = "sequential", status="old", action="read", iostat=iostat)
    
    if(iostat .eq. 0) then
       read(unit = unit, iostat=iostat) data
    else
       print *,"Error, iostat = ", iostat
    end if

    close(unit)
  end subroutine mtxReadBinary_2D

  subroutine mtxReadText_2D(data, Dim, unit, fmt, file, iostat)

    character(*) :: file,fmt
    integer      :: unit, iostat, Dim(2)
    real         :: data(:,:)
    
    open(unit, file = trim(file), form="formatted", access = "sequential", status="old", action="read", iostat=iostat)

    if(iostat .eq. 0) then
       do i=1, Dim(1)
          do j=1, Dim(2)
             read(unit, fmt) data(i,j)
          end do
          read(unit,*) 
       end do
    else
       print *,"Error, iostat = ", iostat
    end if

    close(unit)
  end subroutine mtxReadText_2D


  subroutine mtxWriteBinary_2D(data, unit, file, iostat)
    implicit none

    character(*)   :: file
    real           :: data(:,:)
    integer        :: unit
  
    integer        :: iostat
    
    iostat = 0
    open(unit, file = trim(file), form = "unformatted", access = "sequential", status="replace", action="write", iostat=iostat)

    if(iostat .eq. 0) then 
       write(unit, iostat = iostat) data
    else
       print *,"Error, iostat = ", iostat
    end if
    close(unit)


  end subroutine mtxWriteBinary_2D


 subroutine mtxWriteText_2D(data, Dim, unit, fmt, file, iostat)

    character(*) :: file, fmt
    integer      :: unit, iostat, Dim(2)
    real         :: data(:,:)
    
    open(unit, file = trim(file), form="formatted", access = "sequential", status="replace", action="write", iostat=iostat)
    if(iostat .eq. 0) then
       do i=1, Dim(1)
          do j=1, Dim(2)
             write(unit, fmt) data(i,j)
          end do
          write(unit,*) ""
       end do

    else
       print *,"Error, iostat = ", iostat
    end if

    close(unit)

  end subroutine mtxWriteText_2D



  subroutine mtxReadText_1D(data, dmn, unit, fmt, file, iostat)

    character(*) :: file,fmt
    integer      :: i, unit, iostat, dmn
    real         :: data(:)
    
    open(unit, file = trim(file), form="formatted", access = "sequential", status="old", action="read", iostat=iostat)
    if(iostat .eq. 0) then
       do i=1, dmn
          read(unit, fmt) data(i)
       end do
    else
       print *,"Error, iostat = ", iostat
    end if

    close(unit)
  end subroutine mtxReadText_1D



 subroutine mtxWriteText_1D(data, dmn, unit, fmt, file, iostat)

    character(*) :: file, fmt
    integer      :: i, unit, iostat, dmn
    real         :: data(:)
    
    open(unit, file = trim(file), form="formatted", access = "sequential", status="replace", action="write", iostat=iostat)
    if(iostat .eq. 0) then
       do i=1, dmn
          write(unit, fmt) data(i)
       end do
    else
       print *,"Error, iostat = ", iostat
    end if

    close(unit)

  end subroutine mtxWriteText_1D

end module matrix_tools


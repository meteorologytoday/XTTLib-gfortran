module gnuplot

contains


subroutine mtxWriteGNUPlotBinary_2D(X1, X2, data, unit, file, iostat)

  implicit none
  real           :: X1(:), X2(:), data(:,:)
  character(*)   :: file
  integer        :: unit
  
  integer        :: iostat, i, j, size1, size2
  
  iostat = 0
  open(unit, file = trim(file), form = "unformatted", access = "sequential", status="replace", action="write",&
       iostat=iostat, recl=4)
  
  if(iostat .eq. 0) then 
     size1 = size(X1, dim=1)
     size2 = size(X2, dim=1)

     print *, size1, "*", size2

     write(unit, iostat = iostat) size2


     do j = 1, size2
        write(unit, iostat = iostat) X2(j)
     end do


     write(unit, iostat = iostat) 

     do i = 1, size1
        write(unit, iostat = iostat) X1(i)
        do j = 1, size2
           write(unit, iostat = iostat) data(i,j)
        end do
        write(unit, iostat = iostat) data(i,j)

     end do

  else
     print *,"Error, iostat = ", iostat
  end if

  close(unit)
  
end subroutine mtxWriteGNUPlotBinary_2D


subroutine mtxWriteGNUPlotText_2D(data, X, Y,  Dim, unit, fmt, file, iostat)

  character(*) :: file, fmt
  integer      :: unit, iostat, Dim(2)
  real         :: data(Dim(1), Dim(2)), X(Dim(1), Dim(2)), Y(Dim(1), Dim(2))
  
  open(unit, file = trim(file), form="formatted", access = "sequential", status="replace", action="write", iostat=iostat)
  if(iostat .eq. 0) then
     do i=1, Dim(1)
        do j=1, Dim(2)
           !print *, "(i,j) = (", i, ",",j,"), data = ", data(i,j)
           write(unit, fmt) X(i,j), Y(i,j) ,data(i,j)
        end do
        write(unit,*) ""
     end do
       
  else
     print *,"Error, iostat = ", iostat
  end if
    
  close(unit)
end subroutine mtxWriteGNUPlotText_2D


end module gnuplot

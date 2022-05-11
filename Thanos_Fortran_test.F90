 PROGRAM Thanos_fortran
   IMPLICIT NONE
    real:: y(3,3) = transpose(reshape((/1,2,3,4,5,6,7,8,9/), (/3,3/))),mean(3,3),prueba(10,10),r,vector(9)
    integer:: i,j
    
    
    
    write(*,*) 'matriz'
    write(*,*) y(1,:)
    write(*,*) y(2,:)
    write(*,*) y(3,:)
    
    write(*,*) 'Vector'
    vector=reshape(transpose(y),(/9/),order=(/1/))
    
    write(*,*) vectoccr
!    write(*,*) 'Y====='
!    do i=0,2
!        do j=0,2
!                call getmeanneigh(y,i,j,1,mean(i+1,j+1))
!        end do
!        write(*,*) y(i+1,:)
!    end do
!
!   mean=10+2*(mean-5)/3
!   mean=exp(mean)
!   write(*,*) 'Mean =====' 
!   do i=1,3
!        write(*,*) mean(i,:)
!   end do
!   
!
!   
!   
!   call random_number(r)
!   i=int(r*100)
!   write(*,*) 'Numero 1',i
!   
!      call random_number(r)
!   i=int(r*100)
!   write(*,*) 'Numero 2',i
!   
!      call random_number(r)
!   i=int(r*100)
!   write(*,*) 'Numero 3',i
!   
!      call random_number(r)
!   i=int(r*100)
!   write(*,*) 'Numero 4',i
   
   contains
    SUBROUTINE getmeanneigh(matrix,i,j,r,mean)
     real, INTENT(IN OUT) :: matrix(0:,0:)
     real, INTENT(IN OUT) :: mean
     integer,intent(in):: i,j,r
     integer::  shape_matrix(2),number_neigh
     real :: neighbours
     
     shape_matrix=shape(matrix)
     
     neighbours=sum( matrix( max(0,i-r):min(i+r,shape_matrix(1)-1) , max(0,j-r):min( j+r,shape_matrix(2)-1)  )  )-matrix(i,j)
     number_neigh=size( matrix(max(0,i-r):min(i+r,shape_matrix(1)-1), max(0,j-r):min(j+r,shape_matrix(2)-1)  ) )   -1
     mean=neighbours/number_neigh    
     END SUBROUTINE getmeanneigh
   
  
   END PROGRAM thanos_fortran
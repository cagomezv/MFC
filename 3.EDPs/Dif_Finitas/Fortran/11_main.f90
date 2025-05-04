PROGRAM main
  implicit none
  integer::i,j,N
  integer,allocatable::indx(:)
  real(8)::deter
  real(8),allocatable::A(:,:)
  character(len=8)::ch
  
  !lee archivo input
  open(1,file='matriz.inp',status='old',err=12)
  read(1,*,err=12)ch,N
  allocate(A(N,N),indx(N))
  A=0.d0; indx=0
  read(1,*)ch
  do i=1,N
     read(1,*,err=12)(A(i,j),j=1,N)
  end do
  close(1)

  !subrutina de descomposicion LU
  call LUdescomp(A,N,indx,deter)
  
  !resultados
  write(6,*)'matrices LU'
  do i=1,N
     write(6,"(10(F7.4,1x))")(A(i,j),j=1,N)
  end do
  write(6,*)'intercambio de filas' 
  do i=N,1,-1
     if(indx(i).ne.0)&
          write(6,"(2(i3,1x))")i,indx(i)
  end do
  write(6,"(a,F10.6)")'determinante= ',deter

  STOP
12 write(6,*)'ERROR en el archivo matriz.inp'

END PROGRAM main

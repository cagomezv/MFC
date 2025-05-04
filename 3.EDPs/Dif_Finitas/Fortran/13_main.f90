PROGRAM main
  implicit none
  integer::i,j,N,k
  real(8)::deter
  real(8),allocatable::A(:,:)
  character(len=8)::ch
  
  !lee archivo input
  open(1,file='matriz.inp',status='old',err=12)
  read(1,*,err=12)ch,N
  allocate(A(N,N))
  A=0.d0
  read(1,*)ch
  do i=1,N
     read(1,*,err=12)(A(i,j),j=1,N)
  end do
  close(1)

  !subrutina de sistema de ecuaciones
  call inversionMatriz(A,N)
  
  !resultados
  write(6,*)'Matriz inversa:'
  do i=1,N
     write(6,"(10(F7.4,1x))")(A(i,j),j=1,N)
  end do

  STOP
12 write(6,*)'ERROR en el archivo matriz.inp'

END PROGRAM main

PROGRAM main
  implicit none
  integer::i,j,N
  real(8)::deter
  real(8),allocatable::A(:,:),W(:)
  character(len=8)::ch
  
  !lee archivo input
  open(1,file='matriz.inp',status='old',err=12)
  read(1,*,err=12)ch,N
  allocate(A(N,N),W(N))
  A=0.d0; W=0.d0
  read(1,*)ch
  do i=1,N
     read(1,*,err=12)(A(i,j),j=1,N)
  end do
  read(1,*)ch
  do i=1,N
     read(1,*,err=12)W(i)
  end do
  close(1)

  !subrutina de sistema de ecuaciones
  call sistemaEcuaciones(A,W,N)
  
  !resultados
  write(6,*)'Resultados:'
  do i=1,N
     write(6,"('x(',i2,')= ',F10.6)")i,W(i)
  end do

  STOP
12 write(6,*)'ERROR en el archivo matriz.inp'

END PROGRAM main

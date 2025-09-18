PROGRAM main
  implicit none
  integer::i,j,N,k
  real(8)::deter
  real(8),allocatable::A(:,:),d(:),e(:)
  character(len=8)::ch
  
  !lee archivo input
  open(1,file='matriz.inp',status='old',err=12)
  read(1,*,err=12)ch,N
  allocate(A(N,N),d(N),e(N))
  A=0.d0; d=0.d0; e=0.d0
  read(1,*)ch
  do i=1,N
     read(1,*,err=12)(A(i,j),j=1,N)
  end do
  close(1)

  !subrutina de tridiagonalizacion
  call tridiagon(A,N,d,e)

  !resultados
  write(6,*)'elementos de la diagonal:'
  write(6,"(10(F7.4,1x))")(d(i),i=1,N)
  write(6,*)'elementos fuera de la diagonal:'
  write(6,"(10(F7.4,1x))")(e(i),i=1,N)

  STOP
12 write(6,*)'ERROR en el archivo matriz.inp'

END PROGRAM main

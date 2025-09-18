!Resuelve la ecuación estacionaria de Schroedinger para el oscilador
PROGRAM OsciladorCuantico
  implicit none
  integer::i,j,N,Ni,Nf
  real(8),allocatable::A(:,:),d(:),e(:),tem(:)
  real(8)::h,xi,h_2,Rmin,Rmax
  character(len=8)::ch

  Rmin=-10.d0; Rmax=10.d0

  write(6,*)'los 6 primeros autovalores:'

  N=50; Nf=1000 !puntos en la malla
  do while(N.le.Nf)
     
     h=(Rmax-Rmin)/N; h_2=1.d0/h**2
     
     !define las matrices
     allocate(A(N,N),d(N),e(N),tem(N))
     A=0.d0; d=0.d0; e=0.d0   
     do i=1,N
        xi=Rmin+i*h
        d(i)=2.d0*h_2+xi*xi
        e(i)=-h_2  
        A(i,i)=1.d0  !matriz diagonal
     end do

     !subrutina de diagonalizacion
     call diagotri(d,e,N,A,.true.)

     !ordenos los autovalores desde el mínimo
     do i=1,N
        do j=i+1,N
           if(d(j).lt.d(i))then
              xi=d(i); d(i)=d(j); d(j)=xi
              tem(:)=A(:,i); A(:,i)=A(:,j); A(:,j)=tem(:)
           end if
        end do
     end do
     
     !resultados
     write(6,"(i4,6(2x,F9.6))")N,0.5d0*d(1:6)

     if(N.ge.Nf/2)then
        open(unit=1,file="data_oscilador")
        write(1,"('#',i4,6x,6(F10.6,1x))")N/2,0.5d0*d(1:6)
        !norma de las funciones de onda (calculo la integral)
        tem=0.d0
        do i=1,N
           tem(:)=tem(:)+A(i,:)**2
        end do
        tem=h*tem
        !guardo
        do i=1,N
           xi=Rmin+i*h
           write(1,"(7(F10.5,1x))")xi,A(i,1:6)**2/tem(i)
        enddo
     end if

     deallocate(A,d,e,tem)

     N=2*N
     
  end do

  write(6,*)
  write(6,*)'Autovectores en data_oscilador'
  
END PROGRAM OsciladorCuantico

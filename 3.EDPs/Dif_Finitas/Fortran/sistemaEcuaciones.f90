! A es la matriz de coeficientes, a la salida es la descomposicion LU
! N es el tamanio de la matriz.
SUBROUTINE sistemaEcuaciones(A,W,N)
  IMPLICIT NONE
  
  INTEGER,INTENT(IN)::N
  REAL(8),INTENT(INOUT)::A(N,N),W(N)
  INTEGER::i,j,indx(N)
  REAL(8)::yy(N),xx(N),deter,xt
  
  !descomposicion LU
  indx=0
  call LUdescomp(A,N,indx,deter)
  do i=N,1,-1
     if(indx(i).ne.0)then
        xt=W(i); W(i)=W(indx(i)); W(indx(i))=xt
     end if
  end do

  !hacia adelante
  yy(1)=W(1)
  do i=2,N
     yy(i)=W(i)
     do j=1,i-1
        yy(i)=yy(i)-A(i,j)*yy(j)
     end do
  end do

  !hacia atras
  W=0.d0
  W(N)=yy(N)/A(N,N)
  do i=N-1,1,-1
     W(i)=yy(i)
     do j=i+1,N
        W(i)=W(i)-A(i,j)*W(j)
     end do
     W(i)=W(i)/A(i,i)
  end do


END SUBROUTINE sistemaEcuaciones

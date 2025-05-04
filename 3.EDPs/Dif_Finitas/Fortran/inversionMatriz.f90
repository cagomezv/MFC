! A es la matriz que se quiere invertir, a la salida es la!
! matriz inversa. 
! B es la matriz invertida
! N es el tamanio dela matriz.
SUBROUTINE inversionMatriz(A,N)
  IMPLICIT NONE
  
  INTEGER,INTENT(IN)::N
  REAL(8),INTENT(INOUT)::A(N,N)
  INTEGER::i,j,k,indx(N)
  REAL(8)::deter,LL(N,N),UU(N,N),xt(N)
  
  !descomposicion LU
  indx=0
  call LUdescomp(A,N,indx,deter)

  !inversa de la matriz L
  LL=0.d0
  do i=1,N
     LL(i,i)=1.d0
  end do
  do i=2,N
     do j=1,i-1
        LL(i,j)=-A(i,j)
        do k=j+1,i-1
           LL(i,j)=LL(i,j)-A(i,k)*LL(k,j)
        end do
     end do
  end do

  !inversa de la matriz U
  UU=0.d0
  do i=1,N
     UU(i,i)=1.d0/A(i,i)
     do j=i+1,N
        do k=1,j-1
           UU(i,j)=UU(i,j)-A(k,j)*UU(i,k)
        end do
        UU(i,j)=UU(i,j)/A(j,j)
     end do
  end do

  !calculo la matriz inversa
  A=0.d0
  do i=1,N
     do j=1,N
        do k=1,N
           A(i,j)=A(i,j)+UU(i,k)*LL(k,j)
        end do
     end do
  end do
  !intercambio columnas si se intrcambiaron filas en LU
  do j=N,1,-1
     if(indx(j).ne.0)then
        xt(:)=A(:,j); A(:,j)=A(:,indx(j)); A(:,indx(j))=xt(:)
     end if
  end do

END SUBROUTINE inversionMatriz

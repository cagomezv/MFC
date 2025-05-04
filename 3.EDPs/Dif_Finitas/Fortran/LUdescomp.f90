! A es la matriz que se descompone y LU se almacena en A mismo.
! N es el tamanio dela matriz. d es el determinante.
! indx contiene informacion del intercambio de filas; si indx=0 no se
!     intercambio; el # almacenado se intercambio con el numero del 
!     indice,e.g, indx(2)=5 significa que las filas 2 y 5 se intercam-
!     biaron; se debe intercambiar desde el ultimo hacia el primero 
!     para obtener el A original.
! version modificada desde Numerical Recipes.
SUBROUTINE LUdescomp(A,N,indx,d)
  IMPLICIT NONE
  
  INTEGER,INTENT(IN)::N
  INTEGER,INTENT(OUT)::indx(N)
  REAL(8),INTENT(INOUT)::A(N,N)
  INTEGER::i,j,k,imax
  REAL(8)::sum,umax,tem,d
  REAL(8),ALLOCATABLE::vv(:)
  
  allocate(vv(N))
  d=1.d0; indx=0

  !loop en las columnas
  do j=1,N !la fila 1 de U es la misma que A
     !***paso 3 del algoritmo***
     if(j.ne.1)then
        do i=2,j-1
           sum=a(i,j)
           do k=1,i-1
              sum=sum-a(i,k)*a(k,j)
           enddo
           a(i,j)=sum
        enddo
     end if
     !***paso 4 y 5 del algoritmo***
     umax=0.d0
     do i=j,N
        sum=a(i,j)
        if(j.ne.1)then
           do k=1,j-1
              sum=sum-a(i,k)*a(k,j)
           end do
           a(i,j)=sum
        end if
        !busco el pivot mayor en cada columna
        tem=abs(sum)
        if(tem.gt.umax)then
           umax=tem; imax=i
        end if
     end do
     !Intercambio de filas
     if(j.ne.imax)then
        vv(:)=a(imax,:); a(imax,:)=a(j,:); a(j,:)=vv(:)
        !cambio de paridad para el determinante
        d=-d
        tem=indx(imax); indx(imax)=indx(j); indx(j)=tem
        indx(j)=imax
     endif
     if(a(j,j).eq.0.d0)then
        write(6,*)'La matriz es singular: no existe descomposicion LU'
        stop
     end if
     !***completamos el paso 5 del algoritmo***
     do i=j+1,N
        a(i,j)=a(i,j)/a(j,j)
     enddo
     !Calculamos el valor del determinante
     d=d*a(j,j)
  enddo
  deallocate(vv)

END SUBROUTINE LUdescomp

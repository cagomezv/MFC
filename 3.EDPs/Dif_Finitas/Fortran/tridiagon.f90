! Realiza una reduccion de Housholder a uma matriz simetrica A(:,:)
! A(:,:) es reemplazada por los elementos de la matriz de 
! transformacion S
! d(:) son los elementos de diagonales de la matriz tridiagonal
! e(:) son los elementos no diagonales de la matriz tridiagonal.
!      Siempre e(1)=0 ya que los elementos no diagonales estan 
!      desde i=2..N
! version modificada desde Numerical Recipes (tred2).
SUBROUTINE tridiagon(A,N,d,e)
  IMPLICIT NONE
  INTEGER,INTENT(IN)::N
  REAL(8),INTENT(INOUT)::A(N,N)
  REAL(8),INTENT(OUT)::d(N),e(N)
  integer::i,j,k,l
  real(8)::f,g,h,hh,scale

  !comprueba si la matriz A es simetrica
  do i=1,N-1
     do j=i+1,N
        if(abs(A(i,j)-A(j,i)).gt.1.d-15)then
           write(6,*)'tridiagon: la matriz A no es simetrica'
           stop
        end if
     end do
  end do
  !
  do i=N,2,-1
     l=i-1; h=0.d0; scale=0.d0
     if(l.gt.1)then
        scale=sum(abs(A(i,1:l)))
        if(scale.eq.0.d0)then
           e(i)=A(i,l)
        else
           A(i,1:l)=A(i,1:l)/scale
           h=sum(A(i,1:l)**2)
           f=A(i,l)
           g=-sign(sqrt(h),f)
           e(i)=scale*g
           h=h-f*g
           A(i,l)=f-g
           f=0.d0
           do j=1,l
              A(j,i)=A(i,j)/h
              g=0.d0
              do k=1,j
                 g=g+A(j,k)*A(i,k)
              enddo
              do  k=j+1,l
                 g=g+A(k,j)*A(i,k)
              enddo
              e(j)=g/h
              f=f+e(j)*A(i,j)
           enddo
           hh=f/(h+h)
           do j=1,l
              f=A(i,j)
              g=e(j)-hh*f
              e(j)=g
              do k=1,j
                 A(j,k)=A(j,k)-f*e(k)-g*A(i,k)
              enddo
           enddo
        endif
     else
        e(i)=A(i,l)
     endif
     d(i)=h
  enddo
  d(1)=0.d0; e(1)=0.d0
  do i=1,N
     l=i-1
     if (d(i).ne.0.d0) then
        do j=1,l
           g=0.d0
           do k=1,l
              g=g+A(i,k)*A(k,j)
           enddo
           do k=1,l
              A(k,j)=A(k,j)-g*A(k,i)
           enddo
        enddo
     endif
     d(i)=A(i,i)
     A(i,i)=1.d0
     do j=1,l
        A(i,j)=0.d0; A(j,i)=0.d0
     enddo
  enddo
  
END SUBROUTINE tridiagon

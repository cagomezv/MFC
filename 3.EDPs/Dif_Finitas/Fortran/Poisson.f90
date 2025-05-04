!Calcula la ecuación de Poisson en un sistema bidimensional rectangular
PROGRAM Poisson2D
  implicit none
  real(8),parameter::PI=3.141592653589793d0
  integer::Nx,Ny
  real(8),allocatable::phi(:,:),rho(:,:)
  real(8)::Vx1,Vx2,Vy1,Vy2,Q,eps,Lx,Ly,h,area,rho0
  logical,allocatable::conductor(:,:)
  
  !datos del problema
  write(6,*)'tamaño de la placa rectangular (Lx Ly) en metros?'
  read(5,*)Lx,Ly
  if(Lx.le.Ly)then
     write(6,*)'número de puntos en la grilla en la direccion X?'
     read(5,*)Nx
  else
     write(6,*)'número de puntos en la grilla en la direccionYX?'
     read(5,*)Ny
  end if
  write(6,*)'potencial en las fronteras (Vx1 Vx2 Vy1 Vy2) en Voltios'
  read(5,*)Vx1,Vx2,Vy1,Vy2
  write(6,*)'carga (Q) en Coulumbs?'
  read(5,*)Q
  write(6,*)'presición en la convergencia para la relajación'
  read(5,*)eps
  
  !informacion de la grilla y localiza memoria
  if(Lx.le.Ly)then
     h=Lx/Nx; Ny=int(Ly/h)
  else
     h=Ly/Ny; Nx=int(Lx/h)
  end if
  write(6,*)
  write(6,"(a,i4,'  X ',i4)")'puntos en la grilla: ',Nx,Ny
  allocate(phi(0:Nx,0:Ny),rho(0:Nx,0:Ny),conductor(0:Nx,0:Ny))
  phi=0.d0; rho=0.d0
  area=Lx*Ly
  rho0=Q/area
  rho0=rho0*PI*h**2
  
  !inicializa la red 
  call red_inicial(Vx1,Vx2,Vy1,Vy2,Nx,Ny,rho0,conductor,rho,phi)
  
  !calcula la ecuacion de Poisson iterativamente (Gauss-Seidel)
  call Poisson(Nx,Ny,eps,conductor,rho,phi)
  
  !guarda resultados
  call guardar(Nx,Ny,h,phi)

END PROGRAM Poisson2D

!****************************************************************
SUBROUTINE red_inicial(Vx1,Vx2,Vy1,Vy2,Nx,Ny,rho0,conductor,rho,phi)
  implicit none
  integer,INTENT(IN)::Nx,Ny
  real(8),INTENT(IN)::Vx1,Vx2,Vy1,Vy2,rho0
  real(8),dimension(0:Nx,0:Ny),INTENT(OUT)::rho,phi
  logical,dimension(0:Nx,0:Ny),INTENT(OUT)::conductor
  integer::i,j
  
  !inicia la red con valores 0 and .FALSE.
  phi=0.d0
  conductor=.FALSE.
  rho=0.d0
  !ponemos las condiciones de frontera (conductores)
  conductor(:,0)=.TRUE.  !lado inferior del rectangulo
  conductor(:,Ny)=.TRUE. !lado superior del rectangulo
  conductor(0,:)=.TRUE.  !lado izquierdo del rectangulo
  conductor(Nx,:)=.TRUE. !lado derecho del rectangulo
  phi(:,0)=Vx1
  phi(:,Ny)=Vx2
  phi(0,:)=Vy1
  phi(Nx,:)=Vy2
  
  !distribucion uniforme de carga en el interior
  do i=1,Nx-1
     do j=1,Ny-1
        rho(i,j)=rho0
     enddo
  enddo
  
END SUBROUTINE red_inicial

!****************************************************************
SUBROUTINE Poisson(Nx,Ny,eps,conductor,rho,phi)
  implicit none
  integer,INTENT(IN)::Nx,Ny
  real(8),INTENT(IN)::eps
  logical,dimension(0:Nx,0:Ny),INTENT(IN)::conductor
  real(8),dimension(0:Nx,0:Ny),INTENT(INOUT)::rho,phi
  integer::i,j,iconteo
  real(8)::phi_ij,error,dphi
  
  iconteo=0
  do while (.TRUE.)
     error=0.d0
     do i=1,Nx-1
        do j=1,Ny-1
           !cambiamos el potencial solo para no conductores
           if(.NOT.conductor(i,j))then
              phi_ij=0.25d0*(phi(i+1,j)+phi(i-1,j)+phi(i,j+1)+phi(i,j-1))+ &
                   rho(i,j)
              dphi=abs(phi(i,j)-phi_ij)
              if(error.lt.dphi)error=dphi !error maximo
              phi(i,j)=phi_ij
           endif
        enddo
     enddo
     iconteo=iconteo+1
     if(error.lt.eps)exit
  enddo
  write(6,*)iconteo,' error= ',error
  
END SUBROUTINE Poisson

!****************************************************************
SUBROUTINE guardar(Nx,Ny,h,phi)
  implicit none
  integer,INTENT(IN)::Nx,Ny
  real(8),INTENT(IN)::h
  real(8),dimension(0:Nx,0:Ny),INTENT(IN)::phi
  integer::i,j,n
  real(8)::x,y,V0,v0p,vxy,a,b,pi

  !Solucion analitica para el caso Q=0 con condiciones 0 V0 0 0
  V0=200.d0
  pi=3.141592653589793d0
  v0p=4.d0*V0/pi
  a=Nx*h; b=Ny*h
  !!!!!!!!
  
  open(unit=1,file="data_poisson")
  do i=0,Nx
     do j=0,Ny

        !Solucion analitica para el caso Q=0 con condiciones 0 V0 0 0
        x=i*h; y=j*h; vxy=0.d0
        do n=1,50,2
           vxy=vxy+sin(n*pi*x/b)*sinh(n*pi*y/a)/(n*sinh(n*pi*a/b))
        end do
        vxy=vxy*v0
        !!!!!!!!!!!!
        
        write(1,"(2(F10.5,1x),2(F15.7,1x))")i*h,j*h,phi(i,j),vxy
     enddo
     write(1,*)
  enddo

  write(6,*)'resultados en data_poisson'
  
END SUBROUTINE guardar

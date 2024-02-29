module my_milestones

        use dislin 
        use Linear_systems
        use Cauchy_Problem
        use Temporal_Schemes
        use Fourier_Interpolation
        use Chebyshev_interpolation
        use plots 
        use Boundary_value_problems
        use Collocation_methods
        use Initial_Boundary_Value_Problems
        use Utilities
        use API_Example_Lagrange_Interpolation
        
        use Interpolation
        use Lagrange_interpolation
        use Legendre_points
        
        use API_Example_Cauchy_Problem
        use Numerical_Recipes
        implicit none 

       
    contains  
  


subroutine  milestone_examples
integer :: optionmilestone
real :: p1, p2

     call random_seed ()
     call random_number (p1)
     call random_number (p2)
     
     ! Sólo usar el comando de arriba una vez
     ! Apuntar dichos valores y mantenerlos ya el resto del tiempo

      !p1 = 
      !p2 =
     
     !write(*,*) " ================= Your numbers p1 & p2 are: =================== " 
     !write(*,*) " "
     !write(*,*) " Number p1 = ", p1 
     !write(*,*) " Number p2 = ", p2
     !write(*,*) " "
     write(*,*) " =============================================================== " 
     write(*,*) " ==================== Select a Milestone: ====================== "
     write(*,*) " =============================================================== "
     write(*,*) " "
     write(*,*) " 0. Exit/quit "
     write(*,*) " 1. Milestone 1_2A "
     write(*,*) " 2. Milestone 2B " 
     write(*,*) " 3. Milestone 3 "
     write(*,*) " 4. Milestone 4 "
     write(*,*) " 5. Milestone 5 "
     write(*,*) " 6. Milestone 6 "
     write(*,*) " 7. Milestone 7 "
     write(*,*) " "
     
     read(*,*) optionmilestone
     write(*,*) " "
     
     
     select case(optionmilestone)
         
     case(1)
       call Milestone1_2A  
         
     case(2)   
       call Milestone2B
         
     case(3) 
       call Milestone3 
       
     case(4)
       call Milestone4 
      
     case(5) 
       call Milestone5  
      
     case(6) 
       call Milestone6 
       
     case(7) 

       call Milestone7
       
    case default
          write(*,*) " Case not implemented" 
              
     end select 
    
     write(*,*) " " 
     write(*,*) " =============================================================== "
     write(*,*) " "  
     write(*,*) " Press any key to go Home " 
     read(*,*)
              
end subroutine






!***********************************************************
!* Milestone 1 y 2A
!***********************************************************
subroutine Milestone1_2A 

 ! Usamos la N que nos pidan en el enunciado
 ! Usamos una M razonablemente grande para evitar errores adicionales.
 integer, parameter :: N=20, M=10000
 
 real :: xp1, xp2, xi, p1, p2, comodin
 integer :: x_ev1, x_ev2, x_ev3, x_ev4
 real :: x(0:N), f(0:N), g(0:N)     ! N+1 given points
                                    ! M+1 interpolated points
 
 real :: I_N_f(0:N, 0:M), fe(0:2, 0:M), Error_fe(0:2, 0:M) !f funciton
 real :: I_N_g(0:N, 0:M), ge(0:2, 0:M), Error_ge(0:2, 0:M) !g funciton
 
 real :: Lebesgue_N(-1:N, 0:M),  PI_N(0:N, 0:M)    
 real :: xp(0:M), a=-1, b=1, theta(0:N), theta2(0:M), alpha     
 integer :: i  
 
 !! ============================= Números aleatorios =============================
 
   
 
   
 
    write(*,*) " =============================================================== " 
    write(*,*) " ============== You have selected Milestone 1&2A =============== " 
    write(*,*) " =============================================================== " 
    write(*,*) " "
    write(*,*) " ============ Press enter for Milestone 1&2A results =========== " 
    read(*,*) 
    write(*,*) " Number N = ", N
    write(*,*) " Number M = ", M
    write(*,*) " "
    write(*,*) " =============================================================== " 
    write(*,*) " "

 
 
    
 ! ============================= Puntos de interpolación =============================
 
 ! Nos dirán: Interpolación polinómica global a partir de los puntos interpolación siguientes:
 ! Vendrá expresado con una expresión tal que x_j = ...
 ! Ahí necesitamos meter la expresión que sea, es muy común cos(...)
 
 ! Aquí en theta el argumento de dicho cos(...)
 ! Pueden jugar con llamarlos de las diferentes maneras comentadas más abajo
 ! Pueden pedirnos un valor para una theta distinta a la dada para un sólo apartado
    
 !theta  = [ (PI*i/N, i=0, N) ] !Chebyshev - Gauss - Lobato (extremos)
 
 !theta = [ ((2*i+1)*PI/(2*N+2), i=0, N) ] !Chebyshev - Gauss (ceros)
 
 !theta = [((2*PI*(1+p1)*i)/(2*N*(1+p1)+p2), i = 0, N)]
 
 theta = [((2*PI*i)/(2*N+1), i = 0,N)]
 
 ! Puede no ser cos(...) pero bueno, es lo más normal.
 ! Como contra ejemplo está el Hito 1 (distribución de puntos comentada a continuación)
 
 !x  = [ (a + (b-a)*i/N, i=0, N) ] ! N+1 points
 
 x = cos( theta ) 
 
! También puede ser una función a trozos
! do i = 0,N
!  if(i>N/2) then
!  xi = a + (b-a)*i/N
!  else
!  xi = a + ((b-a)*i/N)**2
!  end if
!  x(i) = xi
! end do
! 
!write(*,*) "A trozos: x = ", x
 
 !!  ============================= Expresión xp ============================= 
 
 !theta2  = [ (PI*i/M, i=0, M) ] !Chebyshev - Gauss - Lobato (extremos)
 
 !theta2 = [ ((2*i+1)*PI/(2*M+2), i=0, M) ] !Chebyshev - Gauss (ceros)
 
 !theta2 = [((2*PI*(1+p1)*i)/(2*M*(1+p1)+p2), i = 0, M)]
 
 theta2 = [ (2*i*PI/(2*M+1), i=0, M) ]
 
 
 !xp = [ (a + (b-a)*i/M, i=0, M) ] ! M+1 points (equispaciado)
  xp = cos( theta2 )
  
 
 !!  ============================= Funciones a estudiar =============================
 
 
 !f = sin ( PI * x )
 !g = 1/( 1 + 25*x**2) 
 
  f = exp(-x**4)
  g = tanh(x)
 
 
 !! ============================= Interpolantes =============================
  
 I_N_f = Interpolant(x, f, N, xp)
 I_N_g = Interpolant(x, g, N, xp)
 
 ! No tocar esto
 Lebesgue_N = Lebesgue_functions( x, xp ) 
 PI_N = PI_error_polynomial( x, xp ) 
 
        
 
 !!  ============================= Funciones y derivadas =============================
 
 ! Spoiler, esto no nos hace falta, sólo vamos a necesitar fe(0,:) y ge(0,:)
 ! En ningún examen de este año 2022/23 ha hecho falta.
 ! fe(1-2,:) y ge(1-2,:) para sacar el error de las derivadas con el interpolante derivado
 
 !fe(0,:) =  sin ( PI * xp ) 
  fe(0,:) = exp(-xp**4)
 !fe(1,:) =  PI * cos ( PI * xp )
 !fe(2,:) =  -PI**2 * sin ( PI * xp ) 
 
 !ge(0,:)= 1/( 1 + 25*xp**2) 
  ge(0,:) =  tanh(xp)
 !ge(1,:) =  -50 * xp / (1 + 25*xp**2)**2  
 !ge(2,:) =  -50  / (1 + 25*xp**2)**2  +(50 * xp)**2 / (1 + 25*xp**2)**2  
 

 
 Error_fe = fe - I_N_f(0:2, :)
 Error_ge = ge - I_N_g(0:2, :) 
 
 
 
 !! ============================= Valores xp de estudio =============================
 
 ! xev es la i con la que consigues el xp del enunciado, la buscamos con lo siguiente:
 ! Ejemplo: Valor del error de interpolación de la función f en x = 0,2 + p2
 ! Usaré xp1 = 0.2+p2, saco su sorrespondiente x_ev1 para poder usarlo luego en las diferentes llamadas a las funciones.
 
 
 !xp1 = -0.90625
 !xp1 = -0.99375
 xp1 = -0.75
 xp2 = 0.55
 
  
 
  !x_ev1 = (xp1 - a)*M/(b-a) !Equispaciada
  !x_ev1 = int(((acos(xp1)*M )/PI)) !Chebyshev - Gauss - Lobato (extremos)
  !x_ev1 = int((((acos(xp1)*(2*M+2))/PI)-1)/2) !Chebyshev - Gauss (ceros)
   x_ev1 = int(((acos(xp1)*(2*M+1))/(2*PI)))
   x_ev2 = int(((acos(xp2)*(2*M+1))/(2*PI)))
 
 write(*,*) " ===================== Numerical results ======================= " 
 write(*,*) " " 
   
  write(*,*) " xp1 =  ", xp1 
  write(*,*) " xp2 =  ", xp2 
  
 !write(*,*) " =================== Funciones f y g =================== "
 write(*,*) " "
 !! Valor funcion f en x = xp1
 !write(*,*) " f en xp1 = ", fe(0,x_ev1)
 !! Valor funcion f en x = xp1
 !write(*,*) " g en xp1 = ", ge(0,x_ev1)
 !write(*,*) " "
 !write(*,*) " ========================= Errores ============================= "
 !write(*,*) " "
 !! Error de interolación de f(x) en x = xp1
 !write(*,*) " E_f(x) en xp1 = ", Error_fe(0,x_ev1)
 !! Error de interolación de g(x) en x = xp1
 !write(*,*) " E_g(x) en xp1 = ", Error_ge(0,x_ev1)
 !
 !! Valor máximo del error de interpolación de f(x)
 !write(*,*) " Error_f(x) max =", maxval(Error_fe(0,:))
 !! Valor máximo del error de interpolación de g(x)
 !write(*,*) " Error_g(x) max =", maxval(Error_ge(0,:))
 !
 !! Posición de x donde E_f(x) es máximo
 !write(*,*) " x para E_f(x) max =", maxloc(Error_fe(0,:))
 !! Posición de x donde E_g(x) es mínimo
 !write(*,*) " x para E_g(x) min =", minloc(Error_ge(0,:))
 !write(*,*) " "
 !write(*,*) " =================== Funciones PI y Lebesgue =================== "
 !write(*,*) " "
 !! Valor funcion de Lebesgue Landa_(N) en x = xp1
 !write(*,*) " Lebesgue Landa_(N) en xp2 = ", Lebesgue_N(0,x_ev2)
 !! Valor constante de Lebesgue Landa_(N)
 !write(*,*) " Constante de Lebesgue(x) (max) =", maxval(abs(Lebesgue_N(0,:))) 
 !! Posición x_max del max. de la función de Lebesgue
 !write(*,*) " xmax de Lebesgue(x) max =", maxval(a+(b-a)*maxloc(abs(Lebesgue_N(0,:)))/M)
 !write(*,*) " "
 !! Valor de la función PI_(N+1) en x = xp1
 !write(*,*) " PI(x) en xp1", PI_N(0,x_ev1)
 !! Valor máximo de la función PI_(N+1)
 !write(*,*) " PI(x) Max. = ", maxval((PI_N(0,:)))
 !! Posición x_min del valor max. PI (N+1)
 !write(*,*) " xmax de PI(x) max =", maxval(a +(b-a)*maxloc(abs(PI_N(0,:)))/M)
 ! write(*,*) " "
 !write(*,*) " ==================== Valores interpolantes ==================== "
 !write(*,*) " "
 !! Valor del interpolante de la función f en x = xp1
 !write(*,*) " I_N_f(x) en xp1 = ", I_N_f(0,x_ev1)
 !! Valor del interpolante de la función g en x = xp1
 !write(*,*) " I_N_g(x) en xp1 = ", I_N_g(0,x_ev1)
 !
 !! Valor de la derivada primera del interpolante de la función f en x = xp1
 !write(*,*) " I_N_f'(x) en xp1 = ", I_N_f(1,x_ev1)
 !! Valor de la derivada primera del interpolante de la función g en x = xp2
 !write(*,*) " I_N_g'(x) en xp1 = ", I_N_g(1,x_ev1)
 !
 !!! Error de interpolación de f'(x) en x = xp1
 !!write(*,*) " E_f'(x) en xp1 = ", Error_fe(1,x_ev1)
 !!! Error de interpolación de g'(x) en x = xp1
 !!write(*,*) " E_g'(x) en xp1 = ", Error_ge(1,x_ev1)
 !
 !write(*,*) " "
 !write(*,*) " ==================== Extras ==================== "
 !write(*,*) " "
 !
 ! comodin = Interpolated_value(xp,I_N_f(0,:),xp1,N)  !interpolante en el punto tambien
 ! write(*,*) "Interpolante con Interpolated_value:", comodin
 ! 
 ! comodin = Interpolated_value(xp,I_N_f(1,:),xp1,N)     !primera derivada interpolante
 ! write(*,*) "Interpolante primera derivada con Interpolated_value:", comodin
 !
 ! comodin = Interpolated_value(xp,Lebesgue_N(0,:),xp1,N)  !lebesgue en el punto
 ! write(*,*) "Lebesgue con Interpolated_value:", comodin
 !
 ! comodin = Interpolated_value(xp,PI_N(0,:),xp1,N)    ! pi en el punto
 ! write(*,*) "PI con Interpolated_value:", comodin

 !write(*,*) " 1. Interpolation error [f(x)-In] = ", Error(0,x_ev)
 !write(*,*) " 3. Valor cte de Lebesgue para extremos de Chebyshev / Valor máximo de la función de Lebesgue: ",maxval(abs(Lebesgue_N(0,:)))
 !write(*,*) " 4. Posición x_max del max. de la función de Lebesgue: ", maxval(a + (b-a)*maxloc(abs(Lebesgue_N(0,:)))/M)
 !write(*,*) " 5. Valor max. del error de interpolación [f(x)-In] = ", maxval(Error(0,:))
 !write(*,*) " 7. Valor max. de la funcion PI (N+1) = ", maxval(PI_N(0,:))
 !write(*,*) " 8. Posición x_min del valor max. PI (N+1) = ", maxloc(PI_N(0,:))
 !write(*,*) " 9. Valor de la funcion PI (N+1) en x=0.765 = ", PI_N(0,:)
 !write(*,*) " 10. Valor de la funcion de Lebesgue en x=0.765 = ", Lebesgue_N(0,:)
 
write(*,*) " 1. PI(x) Max. = ", maxval((PI_N(0,:)))
 write(*,*) " "
write(*,*) " 2. xmax de PI(x) max = ", maxval(a +(b-a)*maxloc(abs(PI_N(0,:)))/M)
 write(*,*) " "
write(*,*) " 3. PI(x) en xp1 = ", PI_N(0,x_ev1)
 write(*,*) " "
write(*,*) " 4. Lebesgue Landa_(N) en xp2 = ", Lebesgue_N(0,x_ev2)
 write(*,*) " "
write(*,*) " 5. Constante de Lebesgue(x) (max) = ", maxval(abs(Lebesgue_N(0,:))) 
 write(*,*) " "
write(*,*) " 6. Error_f(x) max = ", maxval(Error_fe(0,:))
 write(*,*) " "
write(*,*) " 7. Error_g(x) max = ", maxval(Error_ge(0,:))

 !write(*,*) " 8. Error_f(x) max (Equispaciada) = ", maxval(Error_fe(0,:))
 
 !write(*,*) " 9. Error_g(x) max (Extremos de Chebyshev) = ", maxval(Error_ge(0,:))
 
 !write(*,*) " 10. Constante de Lebesgue(x) (max) (Extremos de Chebyshev) = ", maxval(abs(Lebesgue_N(0,:))) 
 
 !write(*,*) " "
 !write(*,*) " ======================== Plot results ========================= "
 !write(*,*) " "
 !write(*,*) " Plot results of milestone 1 y 2A " 
 !write(*,*) " Press enter "
 !read(*,*) 
 !call plot1
 
 
contains 


subroutine plot1

  real :: xmin, xmax, ymin, ymax
  xmin = a; xmax = b; ymin = -2; ymax = 2

  call metafl("xwin")
  call page(4000, 4000)
  call scrmod("reverse")
  
  call disini 
  call winfnt("Courier New Bold")
  call titlin( "Interpolated functions", 3); 
  call graf(xmin, xmax, xmin, (xmax-xmin)/10, ymin, ymax, ymin, (ymax-ymin)/10) 
  call color("blue");  call curve( xp, I_N_f(0,:), M+1)
  call color("red");  call curve( xp, I_N_g(0,:), M+1)
  call incmrk(-1);  call marker(21);
  call color("white"); call curve( x, f, N+1)
  call color("white"); call curve( x, g, N+1)
  call color("white"); call height(80);call title 
  call plot_legends( [ "f", "g" ] ) 
  call disfin
  
  call disini 
  call winfnt("Courier New Bold")
  call titlin( "PI Error function", 3); 
  call setscl(PI_N(0:2,:), M+1, "y") 
  call graf(xmin, xmax, xmin, (xmax-xmin)/10, ymin, ymax, ymin, (ymax-ymin)/10) 
  call color("white");call curve( xp, PI_N(0,:), M+1)
  call color("red"); call curve( xp, PI_N(1,:), M+1)
  call color("blue");  call curve( xp, PI_N(2,:), M+1)
  call color("white"); call height(80);  call title 
  call disfin
  
  call disini 
  call winfnt("Courier New Bold")
  call titlin( "Lebesgue function", 3); 
  call axsscl("log","y"); call labels("log","y"); call labdig(-1, "y")
  call setscl(Lebesgue_N(0,:), M+1, "y") 
  call graf(xmin, xmax, xmin, (xmax-xmin)/10, ymin, ymax, ymin, (ymax-ymin)/10 ) 
  call curve( xp, Lebesgue_N(0,:), M+1)
  call color("white"); call height(80);  call title 
  call disfin
  
  call disini 
  call winfnt("Courier New Bold")
  
  call titlin( "Interpolation Error", 3) 
  call setscl(Error_fe, M+1, "y") 
  call graf(xmin, xmax, xmin, (xmax-xmin)/10, ymin, ymax, ymin, (ymax-ymin)/10 ) 
  call curve( xp, Error_fe(0,:), M+1)
  call color("red"); call curve( xp, Error_fe(1,:), M+1)
  call color("blue");  call curve( xp, Error_fe(2,:), M+1)
  call incmrk(-1);  call marker(21);
  call color("white"); call curve( x, f-f, N+1)
  call color("white"); call height(80);  call title 
  
  call disfin
  
 
end subroutine 
end subroutine 
 

!***********************************************************
!* Milestone 2B  
!***********************************************************
subroutine Milestone2B  
 
 integer, parameter :: N=4, M=400 
 real :: x(0:N), f(0:N), g(0:N)    ! N+1 given points
                                   ! M+1 interpolated points
 
 real :: I_N_f(0:N, 0:M), fe(0:2, 0:M), Error_fe(0:2, 0:M), P_N_f(2, 0:M) !f funciton
 real :: I_N_g(0:N, 0:M), ge(0:2, 0:M), Error_ge(0:2, 0:M), P_N_g(2, 0:M) !g funciton
 
 real :: xp(0:M), a=-1, b=1, theta(0:N)
 real :: c1_f(0:N), c1_g(0:N), c2(0:N), PI = 4 * atan(1d0) 
 integer :: i, k  
  
 write(*,*) " =============================================================== " 
 write(*,*) " =============== You have selected Milestone 2B ================ " 
 write(*,*) " =============================================================== " 
 write(*,*) " "
 write(*,*) " Press enter "
 read(*,*) 
 write(*,*) " "
 
 ! ============================= Puntos de interpolación =============================
 
 ! Nos dirán: Interpolación polinómica global a partir de los puntos interpolación siguientes:
 ! Vendrá expresado con una expresión tal que x_j = ...
 ! Ahí necesitamos meter la expresión que sea, es muy común cos(...)
 
 ! Aquí en theta el argumento de dicho cos(...)
 ! Pueden jugar con llamarlos de las diferentes maneras comentadas más abajo
 ! Pueden pedirnos un valor para una theta distinta a la dada para un sólo apartado
    
 theta  = [ (PI*i/N, i=0, N) ] !Chebyshev - Gauss - Lobato (extremos)
 
 !theta = [ ((2*i+1)*PI/(2*N+2), i=0, N) ] !Chebyshev - Gauss (ceros)
 
 !theta = [((2*PI*(1+p1)*i)/(2*N*(1+p1)+p2), i = 0, N)]
 
 !theta = [((2*PI*i)/(2*N+1), i = 0,N)]
 
 ! Puede no ser cos(...) pero bueno, es lo más normal.
 ! Como contra ejemplo está el Hito 1 (distribución de puntos comentada a continuación)
 
 !x  = [ (a + (b-a)*i/N, i=0, N) ] ! N+1 points
 
 x = cos( theta ) 
 
! También puede ser una función a trozos
! do i = 0,N
!  if(i>N/2) then
!  xi = a + (b-a)*i/N
!  else
!  xi = a + ((b-a)*i/N)**2
!  end if
!  x(i) = xi
! end do
! 
!write(*,*) "A trozos: x = ", x
 
 !!  ============================= Expresión xp (NO TOCAR) ============================= 
  
 xp = [ (a + (b-a)*i/M, i=0, M) ] ! M+1 points
 
 
 !!  ============================= Funciones a estudiar =============================
 
 f = 1/( 1 + 25*x**2) 
 !fe= 1/( 1 + 25*xp**2)
 g = sin ( PI * x ) 
 !ge = sin ( PI * xp ) 
 
 
 
 !! ============================= Interpolantes =============================
  
 I_N_f = Interpolant(x, f, N, xp)
 I_N_g = Interpolant(x, g, N, xp)
 
 
 
 !!  ============================= Funciones y derivadas =============================
 
 ! Spoiler, esto no nos hace falta, sólo vamos a necesitar fe(0,:) y ge(0,:)
 ! En ningún examen de este año 2022/23 ha hecho falta.
 
 fe(0,:)  = 1/( 1 + 25*xp**2) 
 !fe(1,:) =  -50 * xp / (1 + 25*xp**2)**2  
 !fe(2,:) =  -50  / (1 + 25*xp**2)**2  +(50 * xp)**2 / (1 + 25*xp**2)**2  
 
 ge(0,:) =  sin ( PI * xp ) 
 !ge(1,:) =  PI * cos ( PI * xp )
 !ge(2,:) =  -PI**2 * sin ( PI * xp ) 
 
 Error_fe = fe - I_N_f(0:2, :)
 Error_ge = ge - I_N_g(0:2, :) 
 
 
 

 
 c1_f = Chebyshev_transform(f) 
 P_N_f(1, :) = Chebyshev_interpolant( c1_f, xp ) 
 
 c1_g = Chebyshev_transform(f) 
 P_N_g(1, :) = Chebyshev_interpolant( c1_g, xp ) 
 
 c2 = 0
 do k=0, N, 2 
     c2(k) = c_hat(k/2) 
 end do 
 
 P_N_f(2, :) = Chebyshev_interpolant( c2, xp ) 
 
  Error_fe(1,:) = fe(0,:) - I_N_f(0,:)
  Error_fe(2,:) = fe(0,:) - P_N_f(2,:)
!  Error(1,:) = fe - P_N(1,:) 

 write(*,*) " "
 write(*,*) " ======================== Plot results ========================= "
 write(*,*) " "
 write(*,*) " Plot results of milestone 2B " 
 write(*,*) " press enter "
 read(*,*) 
 call plot_milestone2B
 
 
contains

!******************************************************************
!  Cehebyshev expansion. Truncated series. 
!       cos( pi x ) = sum c_m T_m(x) = sum c_m cos( m theta ) 
! 
!  Demonstration:  
!
!  exp( i pi x ) = cos( pi x ) + i sin( pi x ) 
!  Taylor expansion 
!  exp( i pi x ) = sum ( (i pi x)**n / n! ( odd plus even) 
!                = sum  ( (-1)**(2n) pi x )**(2n) / (2n)! + i ....
!
!  cos( pi x ) = sum (-1)**(n) ( pi x )**(2n) / (2n)!
!  Since x = cos( theta ) then, 
!
!  cos( pi x ) = sum (-1)**(n) (pi)**(2n)/(2n)! ( cos theta )**(2n) 
! 
!  Moivre  ( cos theta )**n = ( exp(i theta) + exp(-i theta) )**n / 2**n 
!
!   (x + y)**n = sum from k=0 to k=n  (n | k) x**(n-k) y**k 
!
!   ( cos theta )**(2n) = 1/2**(2n) sum ( 2n | k) exp( i theta( 2n-2k) ) 
!                       = 1/2**(2n-1) sum _0 ^{n} ( 2n | k ) cos ( (2n-2 k) theta ) 
!
! If m = 2n -2 k then, n = m/2 + k 
!
!    c_m = sum from n=m to infinity 2(-1)**n ( PI/2)**(2n) / ( (n-m)! (n+m)! ) 
!
!******************************************************************
function c_hat( m ) result(c) 
  integer, intent(in) :: m
  real :: c 
  
  integer :: n
  real :: f 
  
  c = 0 
  do n=m, 100 
              f =  gamma( real(n-m+1) ) * gamma( real(n+m+1) )
              c = c + 2 * (-1)**n * (PI/2)**(2*n) / f 
  end do 
  if (m==0) c = c/2
  
end function 


subroutine plot_milestone2B 

  real :: xmin, xmax, ymin, ymax
  xmin = a; xmax = b; ymin = -2; ymax = 2

  call metafl("xwin")
  call page(4000, 4000)
  call scrmod("reverse")
  
  call disini 
  call graf(xmin, xmax, xmin, (xmax-xmin)/10, ymin, ymax, ymin, (ymax-ymin)/10) 
  call color("white");  call curve( xp, I_N_f(0,:), M+1)
  call color("blue");   call curve( xp, P_N_f(1,:), M+1)
 ! call color("red");  call curve( xp, P_N(2,:), M+1)
 ! call color("orange");  call curve( xp, fe, M+1)
  
  call incmrk(-1);  call marker(21);
  call color("red"); call curve( x, f, N+1) 
  
  call plot_legends( [ "I_N", "P_N",  "cos(PI x)" ] ) 
  call plot_title( ["Comparison (N=4): Interpolation versus Truncated series", & 
                     "I_N: Chebyshev interpolant, P_N: Truncated series"] ) 
  call disfin
   
  call disini 
  call setscl(Error_fe, M+1, "y") 
  call graf(xmin, xmax, xmin, (xmax-xmin)/10, ymin, ymax, ymin, (ymax-ymin)/10) 
  call color("blue");  call curve( xp, Error_fe(1,:), M+1)
  call color("red");   call curve( xp, Error_fe(2,:), M+1)
  
  call incmrk(-1);  call marker(21);
  call color("white"); call curve( x, f-f, N+1) 
  
  call plot_legends( ["E_I", "E_P", "Points"] ) 
  call plot_title( ["Error comparison: Interpolation versus Truncated series with  N = 4", & 
                    "E_I: Interpolation error, E_P: Truncated series error  " ] )
  call disfin
  
  
end subroutine 

end subroutine 


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

                         !!!!!!!! !!!!!!!! !!!!!!!!!
                         !!    !! !!          !!
                         !!!!!!!! !!!!!       !!
                         !!       !!          !!
                         !!       !!!!!!!! !!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!***********************************************************
!* Milestone 3
!***********************************************************
subroutine Milestone3
 
integer, parameter :: M = 200
integer, parameter :: N1 = 2, N2 = 4 !N1 y N2 como q1 y q2 (recordemos que es para q+1 puntos)
integer :: j, i
real :: x1(0:N1), x2(0:N2), xp(0:M)
real :: a = -1, b = 1
real :: k

real :: chosen_xp
real :: tol 



real, allocatable :: x(:), y(:), PI_N(:, :), Lebesgue_N(:,:) 
! real :: PI = 4 * atan(1.) 
 integer :: N, plot_option 
         

real :: Leb3(-1:N1, 0:M), Leb5(-1:N2, 0:M)
real :: PI3(0:N1,0:M), PI5(0:N2,0:M)
         
 write(*,*) " =============================================================== " 
 write(*,*) " =============== You have selected Milestone 3 ================= " 
 write(*,*) " =============================================================== " 
 write(*,*) " "        
         
do j = 0, N1
             
     k = j
     x1(j) = a + (b - a)*k/N1
            
end do
         
do j = 0, N2
             
     k = j
     x2(j) = a + (b - a)*k/N2
            
end do

xp = [ (a + (b-a)*i/M, i=0, M) ] ! M+1 points
           
! Funciones lebesgue
           
Leb3 = Lebesgue_functions(x1,xp)
Leb5 = Lebesgue_functions(x2,xp)
        
! Funciones de error PI
        
PI3 = PI_error_polynomial(x1,xp)
PI5 = PI_error_polynomial(x2,xp)


        
write(*,*) " ============== PI' with q= 2 ============== "
write(*,*) " "
do j = 0, M
     write(*,*) xp(j), PI3 (1, j)
end do

write(*,*) " "

write(*,*) " ============== PI'' with q= 2 ============== "
write(*,*) " "
do j = 0, M
     write(*,*) xp(j), PI3 (2, j)
end do

write(*,*) " "

write(*,*) " ============== PI' with q= 4 ============== "
write(*,*) " "
do j = 0, M
     write(*,*) xp(j), PI5 (1, j)
end do

write(*,*) " "

write(*,*) " ============== PI'' with q= 4 ============== "
write(*,*) " "
do j = 0, M
     write(*,*) xp(j), PI5 (2, j)
end do

write(*,*) " "

write(*,*) " ============== Lebesgue' with q= 2 ============== "
write(*,*) " "
do j = 0, M
     write(*,*) xp(j), Leb3 (1, j)
end do

write(*,*) " "

write(*,*) " ============== Lebesgue'' with q= 2 ============== "
write(*,*) " "
do j = 0, M
     write(*,*) xp(j), Leb3 (2, j)
end do

write(*,*) " "

write(*,*) " ============== Lebesgue' with q= 4 ============== "
write(*,*) " "
do j = 0, M
     write(*,*) xp(j), Leb5 (1, j)
end do

write(*,*) " "

write(*,*) " ============== Lebesgue'' with q= 4 ============== "
write(*,*) " "
do j = 0, M
     write(*,*) xp(j), Leb5 (2, j)
end do

! Esto por si nos piden para un xp en concreto y no queremos buscarlo a mano, podemos meter una M enorme y modificar la función elegida:

! PI' with q = 2 ---> xp(j), PI3 (1, j) 
! PI'' with q = 2 ---> xp(j), PI3 (2, j)
! PI' with q = 4 ---> xp(j), PI5 (1, j) 
! PI'' with q = 4 ---> xp(j), PI5 (2, j)

! Leb' with q = 2 ---> xp(j), Leb3 (1, j) 
! Leb'' with q = 2 ---> xp(j), Leb3 (2, j)
! Leb' with q = 4 ---> xp(j), Leb5 (1, j) 
! Leb'' with q = 4 ---> xp(j), Leb5 (2, j)

! No olvidar meter una M grande xd (si no, no habrá un punto xp dentro de la tolerancia, aunque podemos cambiarla como es lógico)
! Es una manera cutre de hacerlo de momento que vea si hay alguna función que lo haga o lo piense de otra manera

chosen_xp = 0.982362
tol = 1.0e-7  

write(*,*) " "
write(*,*) " ============== PI' with q= 2 for chosen xp ============== "
write(*,*) " "
write(*,*) " xp =  ", chosen_xp
write(*,*) " "

do j = 0, M
     if (abs(xp(j) - chosen_xp) < tol) then
          write(*,*) xp(j), PI3 (1, j)
     end if
end do


write(*,*) " ================== Plot? 1.Yes 2.No ==================== "
write(*,*) " "

read(*,*) plot_option

select case(plot_option)

case(1)

do N=2, 8, 2 
      
   call dislin_ini(xmax = b, xmin= a, ymax = 15., ymin = 0.0) 
   
   call plot_title( ["Finite difference formulas. First and second derivative" ] )

   allocate( x(0:N), y(0:N), PI_N(0:N,0:M) ) 
   x  = [ (a + (b-a)*i/N, i=0, N) ] 
   y = 0 
 
   PI_N = PI_error_polynomial( x, xp ) 
   
   call color("red");  call marker(0); call incmrk(0) 
   call curve(xp, PI_N(1,:), M+1 )
   
   call color("blue"); call marker(0); call incmrk(0) 
   call curve(xp, PI_N(2,:), M+1 )
   
   
   call color("white");  call marker(21);  call incmrk(-1) 
   call curve(x, y, N+1 ) 
   
   deallocate(x, y, PI_N ) 
   call plot_legends( [ "d PI/dx", "d2 PI/dx",  "points" ] ) 
   call disfin
   
end do 

do N=2, 4, 2 
      
   call dislin_ini(xmax = b, xmin= a, ymax = 30., ymin = 0.0) 
   
   call plot_title( ["Finite difference formulas. First and second derivative" ] )

   allocate( x(0:N), y(0:N), Lebesgue_N(-1:N,0:M) ) 
   x  = [ (a + (b-a)*i/N, i=0, N) ] 
   y = 1 
 
   
   Lebesgue_N = Lebesgue_functions( x, xp ) 
   
   call color("white");  call marker(0); call incmrk(0) 
   call curve(xp, Lebesgue_N(0,:), M+1 ) 
   
   call color("red");  call marker(0); call incmrk(0) 
   call curve(xp, Lebesgue_N(1,:), M+1 )
   
   call color("blue"); call marker(0); call incmrk(0) 
   call curve(xp, Lebesgue_N(2,:), M+1 )
   
   
   call color("white");  call marker(21);  call incmrk(-1) 
   call curve(x, y, N+1 ) 
   
   deallocate(x, y, Lebesgue_N ) 
   call plot_legends( [ "Lebesgue_N", "d Lebesgue_N/dx", "d2 Lebesgue_N/dx",  "points" ] ) 
   call disfin
   
end do 

case(2) 
    end select

end subroutine 




!********************************************************************
!* Milestone4
!*****************************************************************
subroutine Milestone4
 

integer, parameter :: N = 5, M = 1000 , Order = 2
! Cambiar N para la variación de x (Ax), si [0,1] y quieres Ax = 0.1, N = 10 // q te lo dan o te dan el nº de puntos = q+1
real :: x(0:N), xp(0:M), u(0:N), uxk(0:N,2), ErrorUxk(0:N, 2), error_pi(0:N,0:M), Lebesgue_N(-1:N, 0:M), uxk_analytic(0:N, 2), uplot(0:N,2)
real :: x0 = 0, xf = 1 !Intervalo de estudio, hay cambiarlo al que digan
integer :: i
real :: pi = 4 * atan(1.0)
real:: p1, p2


p1 = 0.722973272880877       
p2 = 6.447827978090080E-002

write(*,*) " "
write(*,*) " p1 = ", p1
write(*,*) " p2 = ", p2
write(*,*) " "

 write(*,*) " =============================================================== " 
 write(*,*) " =============== You have selected Milestone 4 ================= " 
 write(*,*) " =============================================================== " 
 
 
    x = [ (x0 + (xf-x0)*i/N, i=0, N) ]
    
    xp = [ (x0 + (xf-x0)*i/M, i=0, M) ]
    
    !!! MUCHISIMO OJO A ESTO NON VS UNIFORM GRID !!!
    call Grid_Initialization( "uniform", "x", x, Order ) ! uso uniform en el examen
    !call Grid_Initialization( "nonuniform", "x", x, Order )
    
    u = EXP(x) / (p1*x**2 + 1.0)
    !u = cosh(x) ! Funcion
    !u = cos(PI*x)
    
     uxk_analytic(0:N, 1) = (p1*x**2 - 2.0*p1*x + 1.0) * EXP(x) / ((p1*x**2 + 1.0)**2) ! Primera derivada
     uxk_analytic(0:N, 2) = (p1**2*x**4 - 4.0*p1**2*x**3 + (6.0*p1**2 + 2.0*p1)*x**2 - 4.0*p1*x - 2.0*p1 + 1.0) * EXP(x) / ((p1*x**2 + 1.0)**3)
    
    
    !uxk_analytic(0:N, 1) = sinh(x) ! Primera derivada
    !uxk_analytic(0:N, 2) = cosh(x) ! Segunda derivada
    !uxk_analytic(0:N, 1) = -PI*sin(PI*x) ! Primera derivada
    !uxk_analytic(0:N, 2) = -PI*PI*cos(PI*x) ! Segunda derivada
    
    !!!!!!!!!!!!!! NORMAL !!!!!!!!!!!!!!!!!!!
    call Derivative( 'x' , 1 , u , uxk(:,1) )
    call Derivative( 'x' , 2 , u , uxk(:,2) )
    ErrorUxk (:,1) = uxk_analytic(:, 1) - uxk(:,1) ! Ojo con el orden para el signo
    ErrorUxk (:,2) = uxk_analytic(:, 2) - uxk(:,2) ! Ojo con el orden para el signo
    
    !!!!! lebesgue and pi errors
    Lebesgue_N = Lebesgue_functions( x, xp ); error_pi = PI_error_polynomial( x,xp )
    !write(*,*) x
    write (*, *) 'Finite differences formulas: ',Order,'th order '
    write (*, '(A, F6.4)') ' Espaciado: ',x(2)-x(1)
    
      !!!! DERIVADAS SEGUNDA !!!!!
    write(*,*) " =================== DERIVADAS SEGUNDAS =================== "
    write(*,*) " "
    write(*,*) "Error 2a derivada f en", x(N/2),"=", ErrorUxk (N/2,2)
    write(*,*) "Error 2a derivada f en", x(N),"=", ErrorUxk (N,2)
    write(*,*) " "
    
    !!!! DERIVADAS PRIMERAS !!!!!
    write(*,*) "=================== DERIVADAS PRIMERAS =================== "
    write(*,*) " "
    write(*,*) "Error 1a derivada f en", x(N/2),"=", ErrorUxk (N/2,1) !Aqui cambias para valores que te pidan (N/2,1), (N/10,1),...
    write(*,*) "Error 1a derivada f en", x(N),"=", ErrorUxk (N,1)
  
 
 
end subroutine 




!********************************************************************
!* Milestone5
!*****************************************************************
subroutine Milestone5


    integer, parameter :: N = 100  ! Numero de puntos = (N+1)
    integer, parameter :: p = 2   ! p = q = orden
    
    integer, parameter ::  Nv = 1, Np = 3
    real :: x(0:N), U(0:N,Nv,Np)
    real :: x0 = 0 , xf = 1 !Intervalo de estudio
    integer :: i, last
    real :: PI = 4*atan(1.) 
    
    real, parameter :: analytical_solution = 0.7634647823
    
    ! Para analizar el error primero tenias que sacar una solucion numerica muy
    ! precisa (N muy grande y q=6 por ejemplo)
    ! Una vez tenias una solucion muy precisa, la tomabas como si fuera la
    ! solucion analitica, y hacias los errores con ella
    
    
    write(*,*) " =============================================================== " 
    write(*,*) " =============== You have selected Milestone 5 ================= " 
    write(*,*) " =============================================================== " 
    write(*,*) " "
    write(*,*) " Solution of boundary value problems "
    write(*,*) " "
    !write(*,*) " y_xx + exp(-x**2) * y_x + - y = 100 * sin(pi*x) * sin(5*pi*x) " 
    !write(*,*) " y(-1) = 0 "
    !write(*,*) " y_x(+1) = 0 "
    !write(*,*) " "
    write(*,*) " Press any key to continue "
    write(*,*) " "
    read(*,*) 
    write(*,*) " =============================================================== " 
    write(*,*) " "
    
    write(*,*) " N = ", N
    
    
    
    x(0) = x0; x(N) = xf  

    !call Grid_Initialization( grid_spacing = "nonuniform", &
    !                         direction = "x",   q = p, nodes = x )
    
    call Grid_Initialization( grid_spacing = "uniform", &
                             direction = "x", q = p, nodes = x )
!   Legendre solution   
    call Boundary_Value_Problem( x = x,                                & 
                                 Differential_operator = ODES,     & 
                                 Boundary_conditions   = BCs, & 
                                 Solution = U(:,:,Np) )
 
 
! call plot(x, transpose(U(:,1,:)),"Milestone5. BVP with order q"  )  
    
    write(*,*) " =============================================================== " 
    write(*,*) " "
    
do i=0,N
   
    last = SIZE(U, DIM=3)
    write(*,*) " x =", x(i), " U(x) =", U(i,1,last)

enddo
 
write(*,*) " "
write(*,*) " =============================================================== " 
write(*,*) " ========== Otras preguntas que pueden hacer en examenes ======= "
write(*,*) " =============================================================== " 
write(*,*) " "
write(*,'(A, F6.4)') " Distancia entre puntos (Delta_x): ", x(3)-x(2)
write(*,*) " "
write(*,*) " Valores para puntos en concreto en función de N "
write(*,*) " "
write(*,*) " x =", x(N/2), " U(x) =", U(N/2,1,last) ! Realmente es más fácil ver la tabla completa que pararse a pensar que N es
write(*,*) " x =", x(0), " U(x) =", U(0,1,last)
write(*,*) " "
write(*,*) " analytical_solution = ", analytical_solution
write(*,*) "              U(N/2) = ", U(N/2,1,last)
write(*,*) " "
write(*,*) " Error at U(N/2) = ", abs(analytical_solution - U(N/2,1,last))
write(*,*) " "
write(*,*) " Recuerda cambiar el parametro 'analytical_solution' para sacar el error "

! Para analizar el error primero tenias que sacar una solucion numerica muy
! precisa (N muy grande y q=6 por ejemplo)
! Una vez tenias una solucion muy precisa, la tomabas como si fuera la
!solucion analitica, y hacias los errores con ella
 
contains 


!****** Differential operator *********
function ODES(x, y, yx, yxx) result(L)
   real, intent(in) :: x, y(:), yx(:), yxx(:)   
   real :: L(size(y)) 
   real :: PI = 4*atan(1d0)

     
   L = yxx + yx + 10*y
   
    !L =  yxx + exp(-x**2) * yx - y - 100*sin(PI*x) * sin(5*PI*x)
    !L = p1*yxx + yx + sin(10*p2*x)*y-cos(20*p1*x) !Metes la funcíon de tu enunciado poniendo todo a un lado de la ecuación
    !L = yxx + yx +30*y-30*sin(PI*x)
    !L = yxx + yx + 10*y
    !L = y*yxx-(yx)**2-1.0
    !L = yxx + exp(-x**2) * yx - y - 100*sin(pi*x) * sin(5*pi*x)
       
end function 
    
!********* Boundary conditions *********
function BCs(x, y, yx) 
           real, intent(in) :: x, y(:), yx(:)  
           real :: BCs(size(y))
           
         
          ! Aquí metemos nuestras condiciones de controno en xo y xf
           
          !BCs = yx - 0    Cuando u' = 0
                           
          !BCs = y - 1     Cuando u = 1
           
          !BCs = y - 1     Cuando u = 1
                           
          !BCs = yx - p1   Cuando u' = p1
           
          ! El que sea que pidan tanto en xo como en xf
           
           

        if (x==x0) then
                           BCs = yx - 1
                           
                           !BCs = y - 1     Cuando u(x0) = 1
                           
                           !BCs = yx - p1   Cuando u'(x0) = p1
                           
        else if (x==xf ) then
            
                           BCs = y - 0                   
        else 
            write(*,*) " Error BCs x=", x; stop  
        endif            
                 
end function  

end subroutine 




!********************************************************************
!* Milestone6
!*****************************************************************
subroutine Milestone6

   integer, parameter :: N = 100 !Time steps
   integer, parameter :: Nv=2, Np=2
   real :: Time(0:N), U(0:N, 2) !U(0:N, Nv, Np) 
   
   real :: a=10., b=28., c=2.6666666666
   real :: t0 = 0, tf = 1, analytical_solution 
   integer :: i, p, last


   
   write(*,*) " =============================================================== " 
   write(*,*) " =============== You have selected Milestone 6 ================= " 
   write(*,*) " =============================================================== " 
   write(*,*) " "
   
   Time = [ (t0 + (tf -t0 ) * i / (1d0 * N), i=0, N ) ]
   
   !Condiciones iniciales
   U(0,1) = 6 ! Funcion
   U(0,2) = 0 ! Derivada
   
   analytical_solution = 0.82327323
  
   ! Poner en Scheme el que te pida:
   ! Adams_Bashforth2
   ! Euler
   ! Crank_Nicolson
   ! Runge_Kutta4
   ! Inverse_Euler
     
        call Cauchy_ProblemS( Time_Domain=Time, Differential_operator=ODE, & 
                              Solution = U, Scheme = Runge_Kutta4 )    

      write(*,*) " "
      write(*,*) " ================= RECUERDA CAMBIAR EL SCHEME! ================= "
      
      write(*,*) " Función "
      write(*,*) " "
      do i=0,N
      write(*,*) " t = ", Time(i), " x'(t) =", U(i,1), " N =", i 
      enddo
      write(*,*) " "  
      write(*,*) " Derivada "
      write(*,*) " "
      do i=0,N
      write(*,*) " t = ", Time(i), " x(t) =", U(i,2), " N =", i 
      enddo
      write(*,*) " "  
        
!write(*,*) " =============================================================== " 
!write(*,*) " ========== Otras preguntas que pueden hacer en examenes ======= "
!write(*,*) " =============================================================== " 
!write(*,*) " "
!write(*,'(A, F6.4)') " Distancia entre tiempos (Delta_t): ", Time(4)-Time(3)
!write(*,*) " "
!write(*,*) " Valores para tiempos en concreto en función de N "
!write(*,*) " "
!write(*,*) " t =", Time(N/2), " x(t) =", U(N/2,1) ! Realmente es más fácil ver la tabla completa que pararse a pensar que N es
!!write(*,*) " t =", Time(0), " x(t) =", U(0,1)
!!write(*,*) " t =", Time(81), " x(t) =", U(81,1) ! t =    2.02500000000000       x(t) =   1.27086609773014       N =          81
!write(*,*) " "
!write(*,*) " analytical_solution = ", analytical_solution
!write(*,*) "              x(N/2) = ", U(N/2,1)
!write(*,*) " "
!write(*,*) " Error at t=(N/2) = ", abs(analytical_solution - U(N/2,1))
!write(*,*) " "
!write(*,*) " raiz(U**2+dU**2) =", sqrt(U(N/10,1)**2+U(N/10,2)**2) !Esto era un ej de un exámen
!write(*,*) " "
!write(*,*) " Valores de la derivada para tiempos en concreto en función de N "
!write(*,*) " "
!write(*,*) " t =", Time(0), "x'(t) =", U(0,2)
!write(*,*) " "
!write(*,*) " Recuerda cambiar el parametro 'analytical_solution' para sacar el error "
!


!   call plot_parametrics(Time, U(:,1,:), ["Euler", "RK4"],"t","x", "Milestone6: Cauchy problem with different time schemes")
!   
!  call Stability_regions_Euler_RK4  
!   
contains

function ODE(U, t) result(F)           
     real :: U(:),t
     real :: F(size(U))
     
     real :: x, dxdt
      
     x = U(1); dxdt = U(2);   

     !F = [ dxdt, cos( 2*t ) - sin(x) ]  !Lo mismo que lo de abajo, primer término no se toca, el segundo términio es despejar u_xx
     
     !U(1) son las x,u (Funcion)
     !U(2) son las xpunto,upunto (Derivadas)
     
     F(1) = U(2) !ESTO NUNCA SE TOCA
       
     !F(2) = -0*U(2)-10*U(1)        !Esto pones a un lado uxx y al otro todo lo demás
     !F(2) = cos(2*t) - sin(U(1))
     F(2) = -U(2)-10*U(1) 
     

 end function
 
end subroutine

subroutine Stability_regions_Euler_RK4
   integer, parameter ::N = 50
   integer :: i, j  
   real :: x(0:N), y(0:N), Region(0:N,0:N)
   real :: x0 = -4d0, xf = 1d0, dx
   real :: y0 = -4d0, yf = 4d0, dy
   
   integer, parameter :: N_levels = 9 
   real :: levels(0:N_levels)  
   
   dx = (xf-x0) / N; dy = (yf-y0) / N
   x = [ ( x0 + dx * i  , i=0,N )]
   y = [ ( y0 + dy * j  , j=0,N )]
   
   levels = [ ( j/real(N_levels)  , j=0, N_levels )]
   do j=1, 2  
     if (j==1)      then 
         call Absolute_Stability_Region(Euler, x, y, Region) 
     else if (j==2) then 
         call Absolute_Stability_Region(Runge_Kutta4, x, y, Region) 
     end if 
     
     call plot_contour(x, y, Region, "Re(z)","Im(z)", levels, graph_type ="isolines"   )  
   end do 
   
end subroutine 




!********************************************************************
!* Milestone7
!*****************************************************************
subroutine Milestone7

       integer, parameter :: Nx = 20, Nt = 1000, Nv = 1
       real ::  x(0:Nx)
       real :: Time(0:Nt), U(0:Nt,0:Nx, Nv)  
       
       real ::  x0 = 0, xf = 1, t0 = 0, tf = 0.2 
       integer :: i, j, k, q = 2 
       integer, parameter :: Nl = 5 
       character(len=10) :: legends(0:Nl) 
     
       write(*,*) " =============================================================== " 
       write(*,*) " =============== You have selected Milestone 7 ================= " 
       write(*,*) " =============================================================== " 
       write(*,*) " "
       
     write (*, '(A50)') 'Milestone 7: Time solution of the 1D heat equation'
     write(*,*) " press enter "
     read(*,*) 
     
     Time = [ (t0 + (tf-t0)*i/Nt, i=0, Nt ) ] 
     do i=0, Nl; legends(i) = "t="//float_to_str(time(Nt/Nl*i)) ; enddo 
         
     x(0) = x0; x(Nx) = xf
     
!    Heat equation 1D  
     call Grid_Initialization( "nonuniform", "x", x, q )
       
     U(0, :, 1)  =  exp(-25*(x-0.5)**2 )
     call Initial_Boundary_Value_Problem(                              & 
                       Time_Domain = Time, x_nodes = x,                & 
                       Differential_operator =  Heat_equation1D,       & 
                       Boundary_conditions   =  Heat_BC1D,             & 
                       Solution = U ) 
     
     call plot_parametrics(x, transpose(U(0:Nt:200,:,1)),           & 
                   legends, "$x$", "$u(x,t)$", "Milestone 7: Heat equation") 
     
     call Stability_Heat_equation_1D
     call Error_Heat1D
   

     end subroutine 

     
     
     
!********************************************************************
!* Functions
!*****************************************************************
function Heat_equation1D( x, t, u, ux, uxx) result(F)
          real, intent(in) ::  x, t, u(:), ux(:), uxx(:)
          real :: F( size(u) ) 
            
            F(1) =   uxx(1)
end function 

function Heat_BC1D(x, t, u, ux) result(BC) 
    real, intent(in) :: x, t, u(:), ux(:) 
    real :: BC( size(u) ) 
    
    real ::  x0 = 0, xf = 1 
    
        if (x==x0) then
                            BC(1) = u(1) 
        else if (x==xf) then
                            BC(1) = u(1)  
        else
             write(*,*)  "Error in Heat_BC1D"; stop 
        endif
end function     
     


!******************************************************************
!   Eigenvalues of the heat equation with FD and Chebyshev 
!
!           d2u/dx2 = lambda u    with u(0)=0 and u(1)=0 
!
!    exact eigenfunctions u_k(x) = sin( k pi x ) with 
!          eigenvalues lambda_k = - ( k pi )**2   
!******************************************************************
subroutine Stability_Heat_equation_1D

   integer, parameter :: Nx = 20
   real ::   x(0:Nx), A(0:Nx, 0:Nx), W(0:Nx, 0:Nx) 
   complex :: lambda(Nx, 0:Nx) ! order=2,4,6,.., index eigenvalue=1,2,3....Nx-1
   real :: rlambda(Nx, 0:Nx) 
   
   integer, parameter ::N = 50, N_levels = 9 
   real :: x_R(0:N), x_I(0:N), dx_R = 5./N, dx_I = 8./N , Region(0:N,0:N), levels(0:N_levels) 
          
   real ::  x0 = 0, xf = 1   
   integer :: i, j, k, q  
   
   x(0) = x0; x(Nx) = xf
  
    
  ! Eigenvalues of discrete operator with different orders 
    do q=2,  Nx, 2 
       call Grid_Initialization( "nonuniform", "x", x, q )
       A =  Linear_operator( 1, x, q, Heat_equation1D, Heat_BC1D )
                
      call Eigenvalues_QR( A, lambda(q,:)) 
      rlambda(q,:) = -Real( lambda(q,:) ) 
      call sortr1( rlambda(q,:), Nx+1, "A") 
      !write(*,*) "lambda =", lambda(q,:) 
      !read(*,*) 
      
    end do 
    
    x_R = [ ( -4 + dx_R * i  , i=0,N )]
    x_I = [ ( -4 + dx_I * j  , j=0,N )]    
    call Absolute_Stability_Region(Runge_Kutta4, x_R, x_I, Region) 
    
    call plot_eigenvalues
    call plot_stability_region 
    call plot_interpolation 
  
           
contains 


subroutine plot_eigenvalues

   real :: xmin, xmax, ymin, ymax
   real :: xi(0:Nx), yi(0:Nx) 
   character(len=10) :: col(15) = ["red","green", "cyan", "orange", "white", "blue", "blue", "cyan", "orange", "white", "green", "blue", "cyan", "orange", "white" ]
   real ::  lambda_exact(0:Nx),  lambda_FD2(0:Nx)  
   
   real :: dx, PI = 4*atan(1.)
   
    dx = (xf-x0)/Nx 
   
    xmin = 0; xmax = Nx; ymin = 0; ymax = 5000;  
    call dislin_ini(xmax, xmin, ymax, ymin) 
   
    do k=0, Nx
      lambda_exact(k)  = (k*PI)**2 
      lambda_FD2(k)  = 4 /dx**2 * ( sin(k*PI*dx/2) )**2 
    end do 
       
    xi = [ (k, k=0,Nx) ]
    do j= 2,  6, 2    
       yi = rlambda(j,:)
       call color(col(j)); call incmrk(1);   call marker(21);  
       call curve( xi, yi, Nx+1 )
    end do 
    
    j=1; yi = rlambda(Nx,:) 
    call color(col(j)); call incmrk(1);   call marker(21);  
    call curve( xi, yi, Nx+1 )
         
    j= 5; yi = lambda_exact
    call color(col(j)); call incmrk(1);   call marker(21);  
    call curve( xi, yi, Nx+1 )
    
    call plot_title( ["Eigenvalues of d2u/dx2 = lambda u ", &
                      " u(0)=0, u(1)=0", " N = 20 " ] ) 
    call plot_legends( ["q=2","q=4", "q=6", "q=N(Chebyshev)","Exact eigenvalues" ] )
    
    call disfin 
 

end subroutine 

subroutine plot_stability_region

   real :: xmin=-4, xmax=1, ymin=-4, ymax=4, dt = 0.001 
   real :: xi(0:Nx), yi(0:Nx) 
       
   
    levels = [ ( j/real(N_levels)  , j=0, N_levels )] 
    call dislin_ini(xmax, xmin, ymax, ymin) 
    
    do i=0, N_levels 
       call contur( x_R, size(x_R), x_I, size(x_I), Region, levels(i) ) 
    end do 
    
    xi = -dt * rlambda(2,:)    
    yi = 0 
    call color("red");  call marker(21);  call incmrk(-1) 
    call curve( xi, yi, Nx+1 )
      
    call plot_title( [" Eigenvalues (q=2) multiplied by the time step (dt=0.001) ", &
                      " inside the stability region of a fourth order Runge-Kutta" ] )
    call disfin 
       

end subroutine 
end subroutine 



subroutine plot_interpolation

 integer, parameter :: N=20, M=500 
 real :: x(0:N), f(0:N)    ! N+1 given points
                           ! M+1 interpolated points 
 real :: xp(0:M), fe(0:M), a=0, b=1, y(0:M, 4)    
 integer :: i, j, k, q(3) = [2, 6, N]
 real, allocatable :: I_N(:,:) 
    
 x  = [ (a + (b-a)*i/N, i=0, N) ] ! N+1 points
 xp = [ (a + (b-a)*i/M, i=0, M) ] ! M+1 points
 
 k = 15
 fe  = sin ( k * PI * xp ) 
 
 do i=1, 3
 call Grid_Initialization( "nonuniform", "x", x, q(i)) 
 f = sin ( k * PI * x ) 
 
 allocate( I_N(0:q(i), 0:M) ) 
 I_N = Interpolant(x, f, q(i), xp) 
 y(:, i) = I_N(0,:) 
 deallocate(I_N) 
end do 
 
 y(:, 4) = fe 
 
 call plot1
 call plot_parametrics(xp, y, ["q=2", "q=6", "q=N", "sin(k pi x)" ], "x", "y", title = "Interpolation with N=20 and k=15")
 
 
 contains
 subroutine plot1

  real :: xmin, xmax, ymin, ymax
  xmin = a; xmax = b; ymin = -2; ymax = 2

  call metafl("xwin")
  call page(4000, 4000)
  call scrmod("reverse")
  
  call disini 
  call winfnt("Courier New Bold")
  call titlin( "Interpolation of sin(k PI x): Chebyshev(blue),  q=6(orange), N=20, k=15", 1) 
  call graf(xmin, xmax, xmin, (xmax-xmin)/10, ymin, ymax, ymin, (ymax-ymin)/10) 
  call color("white");  call curve( xp, y(:,4), M+1) 
  call color("blue");  call curve( xp, y(:,3), M+1)
  call color("orange");  call curve( xp, y(:,2), M+1)
  
  call incmrk(-1);  call marker(21);
  call color("red"); call curve( x, f, N+1)
  
 
  call color("white"); call height(20);call title 
  call disfin
  
 
   
end subroutine 
 
end subroutine





subroutine Error_Heat1D

       integer, parameter :: Nx = 40,  Nv = 1
       real ::  x(0:Nx), U(0:Nx, Nv),  F(0:Nx, Nv), t   
       real ::  R(3, 0:Nx, Nv)
       real ::  x0 = 0, xf = 1
       integer :: i, q, qmax 
  
     do i=1, 2 
         
 !    Polynomial Order=2,4,6, ... Nx           
      q = 2*i 
      
 !    Spatial discretization         
      x(0) = x0; x(Nx) = xf; t = 0 
      call Grid_Initialization( "nonuniform", "x", x, q )
      U = Test_U(Nv, x) 
      F = Spatial_discretization( Heat_equation1D, Heat_BC1D, x, U, t ) 
      
 !    Spatial Truncation Error 
      R(i,:,:) = Spatial_Truncation_Error( Nv, Heat_equation1D, Heat_BC1D, x, q, Test_U  )   
      
     end do  
     R(3,:,:) = Test_Uxx(Nv, x) - F
   
!   call plot(x, F(:,1), "Spatial discretization Heat1D F = Uxx")     
    call plot(x, U(:,1), "Test condition Heat1D U ")   
    call plot(x(1:Nx-1), R(:, 1:Nx-1,1), "Spatial truncation error R with q=2,4 (Richardson & Exact)")
      
     
contains     

function Test_U( Nv, x) result(U)
      integer, intent(in) :: Nv 
      real, intent(in) ::  x(0:) 
      real :: U( 0:size(x)-1, Nv ) 
            
   U(:,1)  =   exp(-80*(x-0.5)**2) 
          
end function 

function Test_Uxx( Nv, x) result(Uxx)
      integer, intent(in) :: Nv 
      real, intent(in) ::  x(0:) 
      real :: Uxx( 0:size(x)-1, Nv ) 
            
 ! U(:,1)   =   exp(-80*(x-0.5)**2) 
 ! Ux(:,1)  =  -160* (x-0.5)  * exp(-80*(x-0.5)**2)  
   Uxx(:,1) = -160 * exp(-80*(x-0.5)**2) -160*(x-0.5) * exp(-80*(x-0.5)**2)  *(-160)* (x-0.5)
          
end function 

end subroutine 






end module  


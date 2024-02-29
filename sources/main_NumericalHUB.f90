!**********************************************************
!* Book:  How to learn Applied maths (amazom.es)
!*        researchgate: Juan A Hernandez Ramos (UPM, Spain)
!**********************************************************    
program main_NumericalHUB

       
       use API_Example_Systems_of_Equations
       use API_Example_Lagrange_Interpolation
       use API_Example_Cauchy_Problem
       use API_Example_Finite_Differences
       use API_Example_Boundary_Value_Problem
       use API_Example_Initial_Boundary_Value_Problem
       use API_Example_IBVP_and_BVP
       use Advanced_problems
       use my_examples
       use my_milestones
       
       use API_Example_IBVP_Chebyshev
       use Special_IBVP
       
       use Boundary_layer
       use Navier_Stokes_cavities
       use Burgers 
       
       use API_Example_Chebyshev_Fourier_interpolation
       implicit none 
       
       
   
       call main 
       
contains 
    
subroutine main
integer :: option 
real :: p1,p2
      
option = 1        
do while (option>0) 
     write(*,*) " =============================================================== " 
     write(*,*) " =========== Welcome to NumericalHUB ft. @aeropedrax =========== " 
     write(*,*) " ==================== This is NumericalHUB+ ==================== " 
     write(*,*) " =============================================================== "
     write(*,*) "  "
     write(*,*) "                            DISCLAIMER:                          "
     write(*,*) "     This is a non-profit remake of the original NumericalHUB    "
     write(*,*) "           Original one is avalible in @jahrWork Github          "
     write(*,*) "             https://github.com/jahrWork/NumericalHUB            "
     write(*,*) " "
     write(*,*) " =============================================================== "
     write(*,*) " ====================== Select an option: ====================== " 
     write(*,*) " =============================================================== "
     write(*,*) " "
     write(*,*) " 0. Exit/quit "
     write(*,*) " 1. Systems of equations "
     write(*,*) " 2. Lagrange interpolation " 
     write(*,*) " 3. Finite difference "
     write(*,*) " 4. ODE Cauchy problems "
     write(*,*) " 5. Boundary value problems "
     write(*,*) " 6. Initial-boundary value problems "
     write(*,*) " 7. Mixed problems: IBVP+BVP "
     write(*,*) " 8. High order ODE schemes "
     write(*,*) " 9. Spectral methods and Navier Stokes examples "
     write(*,*) " 10. My examples "
     write(*,*) " 11. My milestones "
     write(*,*) " "
     write(*,*) " =============================================================== " 
     write(*,*) " "
     
     !read(*,*)
     !read(*,*) option

     !option = 11
     write(*,*) " "

     
     select case(option)
         
     case(0)   
        stop
         
     case(1) 
         call Systems_of_Equations_examples
         
     case(2)   
         call Lagrange_Interpolation_examples 
         
     case(3) 
         call Finite_difference_examples
       
     case(4)
         call Cauchy_problem_examples
      
     case(5) 
         call BVP_examples
      
     case(6) 
         call IBVP_examples 
      
     case(7) 
         call Nonlinear_Plate_Vibrations
       
     case(8) 
         call Advanced_Cauchy_problem_examples
         
     case(9) 
         call Spectral_and_advanced_problems 
         
     case(10) 
         call myExamples
         
      case(11) 
         call milestone_examples   
         
    case default
          write(*,*) "Case not implemented" 
              
      end select 
      
     
end do
  
end subroutine  

end program  


  

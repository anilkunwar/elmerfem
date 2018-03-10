    SUBROUTINE getConductivity( model, n, t, Conductivity )
    ! modules needed
    USE DefUtils

    IMPLICIT None

    ! variables in function header
    TYPE(Model_t) :: model
    INTEGER :: n
    !REAL(KIND=dp) :: dummyArgument
     REAL(KIND=dp) :: t

    ! variables needed inside function
    INTEGER :: DIM
    REAL(KIND=dp) :: x, y,z
    REAL(KIND=dp),POINTER ::  Conductivity(:,:) ! this size needs to be consistent with the sif file!
    !REAL(KIND=dp) :: lambda, beta
    REAL,PARAMETER :: lambda=34.2
    REAL,PARAMETER :: beta=10.0
    Logical :: FirstVisited = .TRUE.
     PRINT *,'size',t,SIZE(Conductivity)

    ! remember these variables
    SAVE DIM, FirstVisited
    FirstVisited=.False.
    x = model % Nodes % x(n)
    y = model % Nodes % y(n)
    z = model % Nodes % z(n)

    
    !Conductivity = 0.0D00
    Conductivity(1,1) = lambda
    Conductivity(1,2) = 0.0D00
    Conductivity(1,3) = lambda*(beta*y)/(x*x+y*y+0.00005)
    Conductivity(2,1) = 0.0D00
    Conductivity(2,2) = lambda
    Conductivity(2,3) = -lambda*(beta*x)/(x*x+y*y+0.00005)
    Conductivity(3,1) = lambda*(beta*y)/(x*x+y*y+0.00005)
    Conductivity(3,2) = (-lambda*beta*x)/(x*x+y*y+0.00005)
    Conductivity(3,3) = lambda*(1.0D3+beta*beta/(x*x+y*y+0.00005))

    END SUBROUTINE getConductivity




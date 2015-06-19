!------------------------------------------------------------------------------
! Subroutine for computing the J-integral
! Peter Råback, 2012
!------------------------------------------------------------------------------
   SUBROUTINE Jintegral( Model,Solver,Timestep,TransientSimulation )
!------------------------------------------------------------------------------
     USE DefUtils

     IMPLICIT NONE
!------------------------------------------------------------------------------
 
     TYPE(Model_t), TARGET :: Model
     TYPE(Solver_t) :: Solver
     REAL(KIND=dp) :: Timestep
     LOGICAL :: TransientSimulation
 !------------------------------------------------------------------------------
!    Local variables
!------------------------------------------------------------------------------
     TYPE(Element_t),POINTER :: CurrentElement
     TYPE(ValueList_t), POINTER :: Params, BC
     INTEGER :: t,dim
     LOGICAL :: Stat
     REAL(KIND=dp) :: Jint,Stot


!------------------------------------------------------------------------------
!    Get variables needed for solution
!------------------------------------------------------------------------------

     Params => Solver % Values
 
     dim = CoordinateSystemDimension()
     IF( dim /= 2 ) THEN
       CALL Fatal('Jintegral','Implemented for 2D case only')
     END IF

!------------------------------------------------------------------------------
!    Do some additional initialization, and go for it
!------------------------------------------------------------------------------

     Jint = 0.0d0
     Stot = 0.0d0
     
     CALL Info( 'Jintegral','-------------------------------------', Level=4 )
     CALL Info( 'Jintegral', 'Computing the J-integral', Level=4 )
     CALL Info( 'Jintegral','-------------------------------------', Level=4 )

     DO t=Solver % Mesh % NumberOfBulkElements+1, &
         Solver % Mesh % NumberOfBulkElements + Solver % Mesh % NumberOfBoundaryElements
       
       CurrentElement => Solver % Mesh % Elements(t)
       
       BC => GetBC( CurrentElement )

       ! only compute integral if the flag is active
       IF( GetLogical( BC,'J integral',Stat ) ) THEN      
         CALL ElementalJintegral( CurrentElement, Jint, Stot )
       END IF
     END DO
     
     WRITE(Message,'(a,ES12.3)') 'Jintegral: ',Jint
     CALL Info( 'Jintegral',Message, Level=4 )

     WRITE(Message,'(a,ES12.3)') 'Line length: ',Stot
     CALL Info( 'Jintegral',Message, Level=4 )

    ! add result to list that is automatically saved by SaveScalars solver, if present
     CALL ListAddConstReal(Model % Simulation,'res: Jintegral',Jint)

!------------------------------------------------------------------------------

 CONTAINS

!------------------------------------------------------------------------------
   SUBROUTINE ElementalJintegral( Element, Jint, Stot)
!------------------------------------------------------------------------------

     TYPE(Element_t), POINTER :: Element
     REAL(KIND=dp) :: Jint,Stot
!------------------------------------------------------------------------------
!    Local variables
!------------------------------------------------------------------------------
     INTEGER :: n
     TYPE(Nodes_t) :: Nodes
     INTEGER, POINTER :: NodeIndexes(:)
     REAL(KIND=dp), POINTER :: Basis(:),dBasisdx(:,:),NodalSxx(:),NodalSyy(:),NodalSxy(:), &
         NodalExx(:), NodalEyy(:), NodalExy(:), NodalU1(:), NodalU2(:)
     REAL(KIND=dp) :: detJ, Sxx,Syy,Sxy,Exx,Eyy,Exy,Normal(3),Traction(2),&
         Wenergy,Integrand,U1x,U2x
     INTEGER :: t,N_Integ,istat
     REAL(KIND=dp) :: s,u,v,w
     TYPE(GaussIntegrationPoints_t), TARGET :: IntegStuff
     REAL(KIND=dp), DIMENSION(:), POINTER :: U_Integ,V_Integ,W_Integ,S_Integ
     LOGICAL :: stat,Visited = .FALSE.
     
     SAVE Nodes, NodalSxx, NodalSxy, NodalSyy, NodalExx, NodalEyy, &
         NodalExy, NodalU1, NodalU2, Basis, dBasisDx, Visited
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
!    Allocate some permanent storage, this is done first time only
!------------------------------------------------------------------------------
     IF ( .NOT. Visited ) THEN
       N = Solver % Mesh % MaxElementNodes

       ALLOCATE( NodalSxx( N ), &
           NodalSyy( N ), &
           NodalSxy( N ), &
           NodalExx( N ), &
           NodalEyy( N ), &
           NodalExy( N ), &
           NodalU1( N ), &
           NodalU2( N ), &
           Nodes % x( N ),   &
           Nodes % y( N ),   &
           Nodes % z( N ),   &
           Basis( N ), &
           dBasisdx( N, 3 ), &
           STAT=istat )
 
       IF ( istat /= 0 ) THEN
         CALL Fatal( 'Jintegral', 'Memory allocation error.' )
       END IF

       Visited = .TRUE.
     END IF


      n = Element % TYPE % NumberOfNodes
      NodeIndexes => Element % NodeIndexes       
      Model % CurrentElement => Element
      
      Nodes % x(1:n) = Solver % Mesh % Nodes % x(NodeIndexes)
      Nodes % y(1:n) = Solver % Mesh % Nodes % y(NodeIndexes)
      Nodes % z(1:n) = Solver % Mesh % Nodes % z(NodeIndexes)
      
!------------------------------------------------------------------------------
! Get the field variables at the boundary element 
!------------------------------------------------------------------------------

! if for you the stress and strains are given as 'Stress 1' etc. then
! you should either change this, or update to a more recent Elmer version.
      CALL GetScalarLocalSolution(NodalSxx, 'Stress_xx')
      CALL GetScalarLocalSolution(NodalSyy, 'Stress_yy')
      CALL GetScalarLocalSolution(NodalSxy, 'Stress_xy')

      CALL GetScalarLocalSolution(NodalExx, 'Strain_xx')
      CALL GetScalarLocalSolution(NodalEyy, 'Strain_yy')
      CALL GetScalarLocalSolution(NodalExy, 'Strain_xy')

      CALL GetScalarLocalSolution(NodalU1, 'Displacement 1')
      CALL GetScalarLocalSolution(NodalU2, 'Displacement 2')

!------------------------------------------------------------------------------
!    Integration stuff
!------------------------------------------------------------------------------
     IntegStuff = GaussPoints( element )

     U_Integ => IntegStuff % u
     V_Integ => IntegStuff % v
     W_Integ => IntegStuff % w
     S_Integ => IntegStuff % s
     N_Integ =  IntegStuff % n

!------------------------------------------------------------------------------
!    Numerical integration over gaussian integration points
!------------------------------------------------------------------------------
     DO t=1,N_Integ

       u = U_Integ(t)
       v = V_Integ(t)
       w = W_Integ(t)

!------------------------------------------------------------------------------
!      Basis function values & derivatives at the integration point
!------------------------------------------------------------------------------
       stat = ElementInfo( Element,Nodes,u,v,w,detJ, Basis,dBasisdx)

       ! Normal vector at integration point 
       Normal = NormalVector( Element,Nodes,u,v)

       ! Stress components at integration point 
       Sxx = SUM( Basis(1:n) * NodalSxx(1:n) )
       Syy = SUM( Basis(1:n) * NodalSyy(1:n) )
       Sxy = SUM( Basis(1:n) * NodalSxy(1:n) )
       
       ! Strain components at integration point
       Exx = SUM( Basis(1:n) * NodalExx(1:n) )
       Eyy = SUM( Basis(1:n) * NodalEyy(1:n) )
       Exy = SUM( Basis(1:n) * NodalExy(1:n) )
       
       ! Displacements at the integration point
       U1x = SUM( dBasisdx(1:n,1) * NodalU1 )
       U2x = SUM( dBasisdx(1:n,1) * NodalU2 )
       
       ! Strain energy and traction at the integration point
       Wenergy = ( Sxx*Exx + Syy*Eyy + 2*Sxy*Exy ) / 2
       Traction(1) = Sxx*Normal(1) + Sxy*Normal(2)
       Traction(2) = Sxy*Normal(1) + Syy*Normal(2)
       
       PRINT *,'t',t,sxx,syy,sxy,exx,eyy,exy,u1x,u2x,wenergy,Traction

       ! Assumes that Einstein's notation really used for the 2nd term
       Integrand = Wenergy*Normal(1) - Traction(1)*U1x - Traction(2)*U2x
       
       ! integration weight 
       s = detJ * S_Integ(t)

       Jint = Jint + s * Integrand 

       ! Just for checking that the line integrates as it should
       Stot = Stot + s
     END DO

!------------------------------------------------------------------------------
   END SUBROUTINE ElementalJintegral
!------------------------------------------------------------------------------

 END SUBROUTINE Jintegral
!------------------------------------------------------------------------------

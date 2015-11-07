    SUBROUTINE AdaptiveSolver(Model, Solver, dt, Transient)

       USE DefUtils
       IMPLICIT NONE
       TYPE(Solver_t) :: Solver
       TYPE(Model_t) :: Model
       REAL(KIND=dp) :: dt
       LOGICAL :: Transient
       !----------------------------------------------------------------
       ! Local variables
       !----------------------------------------------------------------
       TYPE(ValueList_t), POINTER :: SolverParams
        TYPE(Mesh_t), POINTER :: Mesh
        TYPE(Variable_t), POINTER :: Var
       REAL(KIND=dp), POINTER :: WrkPntr(:)
       INTEGER :: SaveResults = 1
       REAL(KIND=dp) :: MaxTimeStepSize = 0.0_dp, MinTimeStepSize, CurrentTimeStepSize=0.0_dp, MinSaveInterval
       REAL(KIND=dp) :: Time, Last_Time = 0.0_dp, Last_Save_Time = 0.0_dp, Duration
       REAL(KIND=dp) :: SaveResults_real
       LOGICAL :: Found, Visited = .FALSE.
       LOGICAL :: AdaptiveTimeStepping, Converged, StopSimulation = .FALSE.
       INTEGER :: istat
        TYPE(Solver_t), POINTER :: iSolver
        REAL(KIND=dp), ALLOCATABLE :: xx(:,:), prevxx(:,:,:)
       INTEGER :: i, j, k, n
       INTEGER :: AdaptiveKeepSmallest, KeepTimeStepSize = 0
       INTEGER :: SteadyConverged
       
       SAVE Last_Time
       SAVE CurrentTimeStepSize
       SAVE KeepTimeStepSize
       SAVE Last_Save_Time
       SAVE Visited
       SAVE xx, prevxx
       SAVE StopSimulation

       CALL Info( 'AdaptiveSolver', '*************************************************' )

       !Read in solver parameters
       !-------------------------
       SolverParams => GetSolverParams()
       IF (.NOT. ASSOCIATED(SolverParams)) CALL FATAL('AdaptiveSolver','No Solver section found')
       
       IF ( .NOT. Visited ) THEN
          ! allocate space for previous result
          n = Model % NumberOfSolvers
          j = 0
          k = 0
          DO i=1,n
             iSolver => Model % Solvers(i)
             IF ( ASSOCIATED( iSolver % Variable  % Values ) ) THEN
                IF ( ASSOCIATED( iSolver % Variable % PrevValues ) ) THEN
                   j = MAX( j, SIZE( iSolver % Variable % PrevValues,2 ) )
                END IF
                k = MAX( k, SIZE( iSolver % Variable % Values ) )
             END IF
          END DO
          ALLOCATE( xx(n,k), prevxx( n,k,j ) )
          ! initialize variables, allocate memory
          Mesh => Model % Meshes
          DO WHILE( ASSOCIATED(Mesh) )
             IF ( Mesh % OutputActive ) EXIT
             Mesh => Mesh % Next
          END DO
          CALL SetCurrentMesh( Model, Mesh )
          NULLIFY(WrkPntr)
          ALLOCATE(WrkPntr(1),STAT=istat)
          IF( istat /= 0 ) CALL Fatal('AdaptiveSolver','Memory allocation error')    
          CALL VariableAdd( Model % Variables, Mesh, Solver, 'AdaptiveTimeStepSize', 1, WrkPntr )
          NULLIFY(WrkPntr)
          ALLOCATE(WrkPntr(1),STAT=istat)
          IF( istat /= 0 ) CALL Fatal('AdaptiveSolver','Memory allocation error 3')    
          CALL VariableAdd( Model % Variables, Mesh, Solver, 'SaveResults', 1, WrkPntr )
       END IF
       
       Time = GetTime()
       
       IF ( StopSimulation ) THEN

          ! Set SaveResults variable to zero
          Var => VariableGet( Model % Variables, 'SaveResults' )
          IF(.NOT. ASSOCIATED(Var)) THEN
             CALL FATAL( 'AdaptiveSolver', 'No memory allocated for SaveResults')       
          END IF
          ! It seems impossible to save an Integer, so we convert it to Real
          Var % Values(1) = 0.0_dp
       
          ! End of simulation was reached at previous timestep - now send stop signal
          WRITE(Message, *)  'End of Simulation! Sending stop signal to solver...'
          CALL ListAddConstReal(Model % Simulation,'Exit Condition',1.0_dp)
          
       ELSE
       
          Converged = .TRUE.
          SaveResults = 1
          DO i=1,Model % NumberOfSolvers
             iSolver => Model % Solvers(i)
             SteadyConverged = iSolver % Variable % SteadyConverged
             IF( ( Time > 0 ) .AND. ( SteadyConverged == 0 ) ) THEN
                Converged = .FALSE.
                CALL Info('AdaptiveSolver', 'No convergence for this timestep')
                EXIT
             END IF
          END DO
          
          AdaptiveTimeStepping = ListGetLogical(SolverParams,'Convergence Adaptive Timestepping',Found)
          
          IF ( AdaptiveTimeStepping ) THEN
             MinTimestepSize = ListGetConstReal( SolverParams, 'Convergence Adaptive Min Timestep', Found )
             IF ( .NOT. Found ) MinTimeStepSize = 1e-15
             AdaptiveKeepSmallest = GetInteger(SolverParams, 'Convergence Adaptive Keep Smallest', Found )
             IF ( .NOT. Found ) AdaptiveKeepSmallest = 2

             IF( Converged .OR. CurrentTimeStepSize - MinTimeStepSize < 1e-10 ) THEN
                ! save current result
                DO i=1,Model % NumberOFSolvers
                   iSolver => CurrentModel % Solvers(i)
                   IF ( ASSOCIATED( iSolver % Variable % Values ) ) THEN
                      n = SIZE( iSolver % Variable % Values )
                      xx(i,1:n) = iSolver % Variable % Values
                      IF ( ASSOCIATED( iSolver % Variable % PrevValues ) ) THEN
                         DO j=1,SIZE( iSolver % Variable % PrevValues,2 )
                            prevxx(i,1:n,j) = iSolver % Variable % PrevValues(:,j)
                         END DO
                      END IF
                   END IF
                END DO
             ELSE
                ! - reset to result for last convergence
                DO i=1,Model % NumberOFSolvers
                   iSolver => CurrentModel % Solvers(i)
                   IF ( ASSOCIATED( iSolver % Variable % Values ) ) THEN
                      n = SIZE(iSolver % Variable % Values)
                      iSolver % Variable % Values = xx(i,1:n)
                      IF ( ASSOCIATED( iSolver % Variable % PrevValues ) ) THEN
                         DO j=1,SIZE( iSolver % Variable % PrevValues,2 )
                         iSolver % Variable % PrevValues(:,j) = prevxx(i,1:n,j)
                         END DO
                      END IF
                   END IF
                END DO
                
                ! - reset time to previous value
                Var  => VariableGet( Model % Solver % Mesh % Variables, 'Time' )
                IF ( ASSOCIATED( Var ) )  Var % Values(1)  = Time - CurrentTimeStepSize
                Time = GetTime()
                Write(Message, *) 'Adaptive Timestepping: Time reset to ', Time
                CALL Info('AdaptiveSolver', Message)
                SaveResults = 0
             END IF
          END IF ! AdaptiveTimeStepping

          ! Get duration
          Duration = ListGetConstReal( SolverParams, 'Duration',Found)
          IF ( .NOT. Found ) Duration = 0.0

          ! Get time minimum saving interval
          MinSaveInterval = ListGetConstReal(SolverParams, 'Min Data Saving Interval',Found)
          IF ( .NOT. Found ) MinSaveInterval = 0.0_dp
          
          IF( SaveResults /= 0 ) THEN
             IF( ( Time > 1e-10 ) .AND. ( Time - Last_Save_Time < MinSaveInterval ) ) THEN
                SaveResults = 0
             ELSE
                Last_Save_Time = Time
             END IF
          END IF
                
          ! Check if duration is exceeded since Time was saved the last time
          ! Duration is multiplied by a factor to avoid problems due to rounding errors
           IF ( Time - Last_Time >= Duration - 1e-10 ) THEN
          ! End of simulation after this timestep! Save one last time...
             SaveResults = 1
             StopSimulation = .TRUE.
          END IF
                
          ! Get time step size
          MaxTimeStepSize = ListGetConstReal(SolverParams, 'Adaptive TimeStep Size',Found)
          IF ( .NOT. Found ) MaxTimeStepSize = 0.0
             
          ! Determine timestep size
          IF ( .NOT. Visited ) THEN   
             CurrentTimeStepSize = MaxTimeStepSize
          ELSE IF ( AdaptiveTimeStepping ) THEN
             IF ( Converged ) THEN
                ! check keepsmallest, increase timestep size if desired
                IF ( KeepTimeStepSize > 0 ) THEN
                   KeepTimeStepSize = KeepTimeStepSize - 1
                ELSE
                   ! try to double the timestep until the timestep set by the user is reached
                   CurrentTimeStepSize = CurrentTimeStepSize * 2
                   IF ( CurrentTimeStepSize > MaxTimeStepSize ) THEN
                      CurrentTimeStepSize = MaxTimeStepSize
                   ELSE
                      WRITE( Message, * ), 'Adaptive Timestepping: Time step set to ', CurrentTimeStepSize
                      CALL Info('AdaptiveSolver', Message)
                   END IF
                END IF
             ELSE
                ! set timestep size to half if MinTimestepSize is not reached
                IF( CurrentTimeStepSize - MinTimeStepSize > 1e-10 ) THEN
                   CurrentTimeStepSize = CurrentTimeStepSize / 2
                   IF ( CurrentTimeStepSize < MinTimeStepSize ) CurrentTimeStepSize = MinTimeStepSize
                   WRITE( Message, * ), 'Adaptive Timestepping: Time step set to ', CurrentTimeStepSize
                   CALL Info('AdaptiveSolver', Message)
                ELSE
                   CALL Warn('AdaptiveSolver','Cannot set smaller time step!')
                   CALL Warn('AdaptiveSolver','Minimum time step might be too big to reach convergence.')
                   CALL Warn('AdaptiveSolver','Continuing anyway...')
                END IF
                ! reset KeepTimeStepSize
                KeepTimeStepSize = AdaptiveKeepSmallest
             END IF
          END IF
             
          ! save TimeStepSize
          Var => VariableGet( Model % Variables, 'AdaptiveTimeStepSize' )
          IF(.NOT. ASSOCIATED(Var)) CALL FATAL( 'AdaptiveSolver', 'No memory allocated for AdaptiveTimeStepSize')       
          Var % Values(1) = CurrentTimeStepSize
                   
          ! Save SaveResults variable
          Var => VariableGet( Model % Variables, 'SaveResults' )
          IF(.NOT. ASSOCIATED(Var)) CALL FATAL( 'AdaptiveSolver', 'No memory allocated for SaveResults')       
          ! It seems impossible to save an Integer, so we convert it to Real
          SaveResults_real = REAL(SaveResults)
          Var % Values(1) = SaveResults_real

          IF ( .NOT. StopSimulation ) THEN
             WRITE(Message, *)  ' SaveResults', SaveResults
          ELSE
             WRITE(Message, *)  'This was the last timestep! Preparing end of Simulation...'
          END IF
       END IF


       CALL Info( 'AdaptiveSolver', Message, level=4)
       CALL Info( 'AdaptiveSolver', '*************************************************' )
       
       IF ( .NOT. Visited ) Visited = .TRUE.
       
    END SUBROUTINE AdaptiveSolver

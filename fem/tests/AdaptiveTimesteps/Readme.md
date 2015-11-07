This is the code for Adaptive Timestepping developed by Matthias Zenker (2012).
Additional information can be obtained from 
http://www.elmerfem.org/forum/viewtopic.php?f=3&t=2674&sid=1babe4551462d20d63cee8d04a451f5c
When this solver has to be used the following code has to be put into the Result Output Solver
###############################################################################################
      ! Check SaveResults variable from Adaptive solver, don't save if it is equal to 0
      IF( ListGetLogical(Params,'Check SaveResults Variable',Found) ) THEN
         Var => VariableGet( Model % Variables, 'SaveResults' )
         IF(.NOT. ASSOCIATED(Var)) THEN
            CALL Warn( 'ResultOutputSolver', 'SaveResults not found' )
            SaveResults = 1.0_dp
         ELSE
            SaveResults = Var % Values(1)
         END IF
         IF( SaveResults < 0.1 ) THEN
           CALL Info( 'ResultOutputSolver', 'Nothing to save...')
           CALL Info( 'ResultOutputSolver', '-------------------------------------')
           RETURN
        END IF
      END IF 
#####################################################################################################
and the following has to be written in the solver input file (SIF)
######################################################################################################
    Simulation
    ...
      Timestep intervals = 10000
      Timestep Size = Equals AdaptiveTimeStepSize
    ...
    End
    ...
    Solver 1
      Equation = Adaptive Solver
      Convergence Adaptive Timestepping = Logical True
      Convergence Adaptive Min Timestep = Real 0.1
      Convergence Adaptive Keep Smallest = Integer 5
      Duration = Real 1200
      Variable = -dofs 1 -global AdaptiveTimeStepSize
      Procedure = "AdaptiveSolver" "AdaptiveSolver"
      Adaptive TimeStep Size = 5
      Min Data Saving Interval = Real 1
      Exec Solver = After Timestep
    End

    Solver 2
      Equation = Result Output
      Procedure = "ResultOutputSolve" "ResultOutputSolver"
    ...
      Check SaveResults Variable = Logical True
      Exec Solver = After Timestep
    End
#############################################################################################################

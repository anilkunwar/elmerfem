    !-----------------------------------------------------
    ! material property user defined function for ELMER:
    ! Viscosity of two-phase and liquid eutectic Pb–Bi alloy at high temperatures
    ! S V Stankus 1 , R A Khairulin 1 , and A G Mozgovoy 2
    ! V. Sobolev / Journal of Nuclear Materials 362 (2007) 235–247
    ! visc = eta_0 * exp (E_n/RT), eta_0 = 4.94*10^-4, E_n = 6270.
    ! Mushy Region (Two phase liquid solid) region.
    ! ! asymptotic-visc = 0.091* exp(4.32*f_s)
    ! L.S. Turng, K.K. Wang, J. Mater. Sci. 26 (1991), 2173-2183.
    ! Written By: Anil Kunwar (2017-12-20 Wednesday & 2017-12-21 Thursday)
    !-----------------------------------------------------
    FUNCTION getViscosity( model, n, temp ) RESULT(viscosity)
    ! modules needed
    USE DefUtils
    IMPLICIT None
    ! variables in function header
    TYPE(Model_t) :: model
    INTEGER :: n
    REAL(KIND=dp) :: temp, viscosity
    Real, Parameter::R=8.31 !Universal Gas Constant
    !Real, Parameter::fs=1 !value of solid fraction
    Real, Parameter::amplitude =0.091 !value of ampt term in asymptotic viscosity
   
    !fs is considered as a function of temperature
    ! variables needed inside function
    REAL(KIND=dp) :: eta0, qa, meltT, solidT
    Logical :: GotIt
    TYPE(ValueList_t), POINTER :: material

    ! get pointer on list for material
    material => GetMaterial()
    IF (.NOT. ASSOCIATED(material)) THEN
    CALL Fatal('getViscosity', 'No material found')
    END IF

    ! read in reference viscosity
    eta0 = GetConstReal( material, 'reference viscosity',GotIt)
    IF(.NOT. GotIt) THEN
    CALL Fatal('getViscosity', 'reference viscosity in Viscosity Equation not found')
    END IF

    ! read in activations energy
    qa = GetConstReal( material, 'Activation Energy', GotIt)
    IF(.NOT. GotIt) THEN
    CALL Fatal('getViscosity', 'Activation Energy not found')
    END IF

        ! read in solid fraction
    !fs = GetConstReal( material, 'Solid fraction', GotIt)
    !IF(.NOT. GotIt) THEN
    !CALL Fatal('getViscosity', 'Solid fraction not found')
    !END IF

            ! read in melting temperature
    meltT = GetConstReal( material, 'Tm-rho', GotIt)
    IF(.NOT. GotIt) THEN
    CALL Fatal('getViscosity', 'Melting temperature not found')
    END IF

    
    ! read in solidification temperature
    solidT = GetConstReal( material, 'Ts-rho', GotIt)
    IF(.NOT. GotIt) THEN
    CALL Fatal('getViscosity', 'Solidification temperature not found')
    END IF

    ! compute density conductivity viscosity
    IF (meltT <= temp) THEN ! check for physical reasonable temperature
       CALL Warn('getViscosity', 'The Pb-Bi material is in liquid state.')
    	viscosity = eta0*exp(qa/(R*temp))
    ELSE IF (solidT >= temp) THEN
    	 CALL Warn('getViscosity', 'The Pb-Bi material is in solid state.')
    	viscosity = 6.8421	!PaS
    ELSE 
    	 CALL Warn('getViscosity', 'The Pb-Bi material is in mushy zone.')
    	!viscosity = amplitude*exp(4.32*fs)
        viscosity = amplitude*exp(-2.12112*temp + 844.1064)
     !linear fitting with temperature and fs at mushy zone, fs = -0.491*T + 195.395 for eutectiv Pb-Bi alloy 
    END IF

    END FUNCTION getViscosity


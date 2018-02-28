    !-----------------------------------------------------
    ! material property user defined function for ELMER:
    ! Thermal conductivity of Pb-Bi Eutectic as a function of temperature
    ! (kth_sn)solid = A + B*T + C*T^2, where A = 3.61 and B = 1.517*10^-2, C =-1.741*10^-6
    ! V. Sobolev, J. Nuclear Mater. (2007), Vol. 362:235-247
    ! Written By: Anil Kunwar (2015-03-13)
    !-----------------------------------------------------
    FUNCTION getThermalConductivity( model, n, temp ) RESULT(thcondt)
    ! modules needed
    USE DefUtils
    IMPLICIT None
    ! variables in function header
    TYPE(Model_t) :: model
    INTEGER :: n
    REAL(KIND=dp) :: temp, thcondt

    ! variables needed inside function
    REAL(KIND=dp) :: refThCond, alpha, beta, refTemp
    Logical :: GotIt
    TYPE(ValueList_t), POINTER :: material

    ! get pointer on list for material
    material => GetMaterial()
    IF (.NOT. ASSOCIATED(material)) THEN
    CALL Fatal('getThermalConductivity', 'No material found')
    END IF

    ! read in reference density at reference temperature
    refThCond = GetConstReal( material, 'Constant term',GotIt)
    IF(.NOT. GotIt) THEN
    CALL Fatal('getThermalConductivity', 'Constant term not found')
    END IF

    ! read in term B
    alpha = GetConstReal( material, 'Term B', GotIt)
    IF(.NOT. GotIt) THEN
    CALL Fatal('getThermalConductivity', 'Term B not found')
    END IF
    
    ! read in term C
    beta = GetConstReal( material, 'Term C', GotIt)
    IF(.NOT. GotIt) THEN
    CALL Fatal('getThermalConductivity', 'Term C not found')
    END IF


    ! read in reference temperature
    refTemp = GetConstReal( material, 'Tm-rho', GotIt)
    IF(.NOT. GotIt) THEN
    CALL Fatal('getThermalConductivity', 'Melting temperature not found')
    END IF


    ! compute density conductivity
    IF (refTemp <= temp) THEN ! check for physical reasonable temperature
       CALL Warn('getThermalConductivity', 'The Pb-Bi eutectic material is in liquid state.')
            !CALL Warn('getThermalConductivity', 'Using density reference value')
    !thcondt = 1.11*(refThCond + alpha*(temp))
    thcondt = refThCond + alpha * temp + beta * (temp)**2
    ELSE
    thcondt = refThCond + alpha * temp + beta * (temp)**2
    END IF

    END FUNCTION getThermalConductivity


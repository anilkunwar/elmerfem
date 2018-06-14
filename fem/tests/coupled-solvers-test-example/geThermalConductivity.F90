    !-----------------------------------------------------
    ! material property user defined function for ELMER:
    ! Thermal conductivity of tin as a function of temperature
    ! (kth_sn)solid = A*T + B, where A = -0.06028 and B = 88.9342
    ! (kth_sn)liquid = 1.11*(A*T + B)
    ! Fatma Meydaneri et al. Met. Mater. Int. (2012), Vol. 18:77-85
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
    REAL(KIND=dp) :: refThCond, alpha,refTemp
    Logical :: GotIt
    TYPE(ValueList_t), POINTER :: material

    ! get pointer on list for material
    material => GetMaterial()
    IF (.NOT. ASSOCIATED(material)) THEN
    CALL Fatal('getThermalConductivity', 'No material found')
    END IF

    ! read in reference density at reference temperature
    refThCond = GetConstReal( material, 'Reference Thermal Conductivity B',GotIt)
    IF(.NOT. GotIt) THEN
    CALL Fatal('getThermalConductivity', 'Reference Thermal Conductivity not found')
    END IF

    ! read in Temperature Coefficient of Resistance
    alpha = GetConstReal( material, 'slope A', GotIt)
    IF(.NOT. GotIt) THEN
    CALL Fatal('getThermalConductivity', 'slope of thermal conductivity-temperature curve not found')
    END IF

    ! read in reference temperature
    refTemp = GetConstReal( material, 'Melting Point Temperature', GotIt)
    IF(.NOT. GotIt) THEN
    CALL Fatal('getThermalConductivity', 'Reference Temperature not found')
    END IF


    ! compute density conductivity
    IF (refTemp <= temp) THEN ! check for physical reasonable temperature
       CALL Warn('getThermalConductivity', 'The Sn material is in liquid state.')
            !CALL Warn('getThermalConductivity', 'Using density reference value')
    !thcondt = 1.11*(refThCond + alpha*(temp))
    thcondt = 1.2*(refThCond + alpha*(temp))
    ELSE
    thcondt = refThCond + alpha*(temp)
    END IF

    END FUNCTION getThermalConductivity


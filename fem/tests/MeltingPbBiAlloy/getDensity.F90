    !-----------------------------------------------------
    ! material property user defined function for ELMER:
    ! Density of liquid eutectic Pb–Bi alloy at high temperatures
    ! S V Stankus 1 , R A Khairulin 1 , and A G Mozgovoy 2
    ! ρ(T) = 10524.6 – 1.3571*(T-397.9) + 1.69*10^-4 *(T-397.9)^2
    ! density = A + B* (T-meltT) + c* (T-meltT)^2
    ! refdens = A, term1 = B, term2 = c, meltT = 397.9
    ! T_{mp} = (397.9±0.2) K is the melting point of Pb-Bi eutectic alloy
    ! Written By: Anil Kunwar (2017-12-20 Wednesday & 2017-12-21 Thursday)
    !-----------------------------------------------------
    FUNCTION getDensity( model, n, temp ) RESULT(density)
    ! modules needed
    USE DefUtils
    IMPLICIT None
    ! variables in function header
    TYPE(Model_t) :: model
    INTEGER :: n
    REAL(KIND=dp) :: temp, density

    ! variables needed inside function
    REAL(KIND=dp) :: refdens, term1, term2, meltT
    Logical :: GotIt
    TYPE(ValueList_t), POINTER :: material

    ! get pointer on list for material
    material => GetMaterial()
    IF (.NOT. ASSOCIATED(material)) THEN
    CALL Fatal('getDensity', 'No material found')
    END IF

    ! read in reference density (term 0 or coefficient A) at reference temperature
    refdens = GetConstReal( material, 'CoeffA',GotIt)
    IF(.NOT. GotIt) THEN
    CALL Fatal('getDensity', 'Term A in Density Equation not found')
    END IF

    ! read in term 1 or coefficient B
    term1 = GetConstReal( material, 'CoeffB', GotIt)
    IF(.NOT. GotIt) THEN
    CALL Fatal('getDensity', 'Term B in Density Equation not found')
    END IF

      ! read in term 2 or coefficient c
    term2 = GetConstReal( material, 'Coeffc', GotIt)
    IF(.NOT. GotIt) THEN
    CALL Fatal('getDensity', 'Term c in Density Equation not found')
    END IF


    ! read in the melting temperature
    meltT = GetConstReal( material, 'Tm-rho', GotIt)
    IF(.NOT. GotIt) THEN
    CALL Fatal('getDensity', 'melting temperature in Density Equation not found')
    END IF


    ! compute density conductivity
    IF (meltT <= temp) THEN ! check for physical reasonable temperature
       CALL Warn('getThermalConductivity', 'The Pb-Bi material is in liquid state.')
            !CALL Warn('getThermalConductivity', 'Using density reference value')
    !thcondt = 1.11*(refThCond + alpha*(temp))
    density = refdens - term1 * (temp - meltT) + term2 * (temp-meltT)**2
    ELSE
    density = 10656.1 !Khairulin2005, Journal of Alloys and Compounds, https://doi.org/10.1016/j.jallcom.2004.06.045
    END IF

    END FUNCTION getDensity


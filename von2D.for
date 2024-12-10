C=======================================================================
C   An enhanced local damage model for 2D and 3D quasi-brittle fracture:
C   ABAQUS-FEM implementation and comparative study on the effect 
C   of equivalent strains.
C=======================================================================
SUBROUTINE UMAT(STRESS,STATEV,DDSDDE,SSE,SPD,SCD,
     1 RPL,DDSDDT,DRPLDE,DRPLDT,
     2 STRAN,DSTRAN,TIME,DTIME,TEMP,DTEMP,PREDEF,DPRED,CMNAME,
     3 NDI,NSHR,NTENS,NSTATV,PROPS,NPROPS,COORDS,DROT,PNEWDT,
     4 CELENT,DFGRD0,DFGRD1,NOEL,NPT,LAYER,KSPT,KSTEP,KINC,ELEM_COORDS,AREA)

C-----------------------------------------------------------------------
C   INCLUDE AND VARIABLE DECLARATIONS
C-----------------------------------------------------------------------
    INCLUDE 'ABA_PARAM.INC'

    CHARACTER*80 CMNAME
    DIMENSION STRESS(NTENS), STATEV(NSTATV),
     1 DDSDDE(NTENS,NTENS), DDSDDT(NTENS), DRPLDE(NTENS),
     2 STRAN(NTENS), DSTRAN(NTENS), TIME(2), PREDEF(1), DPRED(1),
     3 PROPS(NPROPS), COORDS(2,4), DROT(3,3), DFGRD0(3,3), DFGRD1(3,3)

    DOUBLE PRECISION CC(NTENS,NTENS), delta(NTENS),
     $ II1, JJ2, II1_d(NTENS), JJ2_d(NTENS), NODEE(1),
     $ AA, e_tilde1, kappa, kappa_n, omega, ee, gg, coeff_omega, coeff_ee(NTENS),
     $ a1, a2, a4(NTENS), Betaa(1), ELEM_COORDS(2,4),
     $ STRESS_EQ(NTENS), CCC(4,4), X(4), Y(4),
     $ dk_dee(1)

C   MATERIAL PROPERTIES
    DOUBLE PRECISION E, NU, k, kappa_o

C   PARAMETERS AND CONSTANTS
    INTEGER i, j, a, b
    PARAMETER(Zero=0d0, One=1.0d0, Two=2.0d0, Three=3.0d0,
     $ Six=6.0d0, TOL=1d-25, Alpha=1.0d0, Beta=305.0d0, Ne=3576.0d0,
     $ Eta=5.0d0)

C-----------------------------------------------------------------------
C   EXTRACT PROPERTIES FROM INPUT
C-----------------------------------------------------------------------
    E = PROPS(1)          ! Young's modulus
    NU = PROPS(2)         ! Poisson's ratio
    k = PROPS(3)          ! Ratio f_t/f_c
    kappa_o = PROPS(4)    ! Damage initiation threshold

C-----------------------------------------------------------------------
C   INITIALIZATION
C-----------------------------------------------------------------------
    CC = 0d0
    delta = 0d0
    DDSDDE = 0d0

C   TOTAL STRAIN INCREMENT
    DO a = 1, NTENS
        STRAN(a) = STRAN(a) + DSTRAN(a)
    END DO

C-----------------------------------------------------------------------
C   ELASTIC STIFFNESS MATRIX AND INVARIANTS (PLANE STRAIN)
C-----------------------------------------------------------------------
    IF (NDI == 3 .AND. NSHR == 1) THEN
        ! Elastic stiffness CC
        DO i = 1, NDI
            DO j = 1, NDI
                CC(i, j) = NU
            END DO
        END DO

        DO i = 1, NDI
            CC(i, i) = 1.0d0 - NU
        END DO

        CC(4, 4) = (One - Two * NU) / Two
        CC = CC * E / ((One + NU) * (One - Two * NU))

        ! Kronecker delta
        delta(1) = One
        delta(2) = One
        delta(3) = One
        delta(4) = Zero

        ! Invariants II1 and JJ2
        II1 = 0d0
        DO a = 1, NTENS
            II1 = II1 + STRAN(a) * delta(a)
        END DO
        JJ2 = STRAN(1)**2 / Three + STRAN(2)**2 / Three - STRAN(1) * STRAN(2) / Three +
     $       (STRAN(4) / Two)**2
    ENDIF

C-----------------------------------------------------------------------
C   EQUIVALENT STRAIN (MODIFIED VON MISES)
C-----------------------------------------------------------------------
    AA = ((((k - One) * II1) / (One - Two * NU))**Two) +
     $    (12.0d0 * k * JJ2) / ((One + NU)**Two)
    ee = ((k - One) * II1) / (Two * k * (One - Two * NU)) +
     $    (SQRT(AA)) / (Two * k)

C   HISTORY VARIABLE kappa_n
    IF (TIME(2) .EQ. 0d0) THEN
        kappa_n = kappa_o
    ELSE
        kappa_n = STATEV(1)
    ENDIF
    STATEV(1) = kappa_n

C   UPDATE HISTORY PARAMETER kappa
    IF (kappa_n .LT. ee) THEN
        kappa = ee
        STATEV(1) = ee
        dk_dee(1) = One
    ELSE
        kappa = kappa_n
        dk_dee(1) = Zero
    ENDIF

C-----------------------------------------------------------------------
C   DAMAGE VARIABLE AND DERIVATIVES
C-----------------------------------------------------------------------
    NODEE = CELENT
    Betaa = E * kappa_o * NODEE / (64.0d0 / 1000.0d0)
    omega = One - (kappa_o / kappa) * ((One - Alpha) +
     $        Alpha * EXP(-Betaa(1) * (kappa - kappa_o)))

    coeff_omega = (kappa_o / (kappa**Two)) * (One - Alpha +
     $        Alpha * EXP(-Betaa(1) * (kappa - kappa_o))) +
     $        (kappa_o / kappa) * (Alpha * EXP(-Betaa(1) * (kappa - kappa_o))
     $        * Betaa(1))

C-----------------------------------------------------------------------
C   STRESS UPDATE (WITH DAMAGE)
C-----------------------------------------------------------------------
    STRESS_EQ = 0d0
    STRESS = 0d0
    DO i = 1, NTENS
        DO j = 1, NTENS
            STRESS_EQ(i) = STRESS_EQ(i) + CC(i, j) * STRAN(j)
            STRESS(i) = (One - omega) * STRESS_EQ(i)
        END DO
    END DO

C-----------------------------------------------------------------------
C   COMPUTE DDSDDE (TANGENT MODULUS)
C-----------------------------------------------------------------------
    CCC = 0d0
    DO i = 1, NTENS
        DO j = 1, NTENS
            CCC(i, j) = STRESS_EQ(i) * a4(j)
        END DO
    END DO

    DO i = 1, NTENS
        DO j = 1, NTENS
            DDSDDE(i, j) = (One - omega) * CC(i, j) -
     $                     coeff_omega * dk_dee(1) * CCC(i, j)
        END DO
    END DO

C-----------------------------------------------------------------------
C   STORE STATE VARIABLES
C-----------------------------------------------------------------------
    STATEV(2) = omega
    STATEV(3) = coeff_omega
    RETURN
END SUBROUTINE UMAT

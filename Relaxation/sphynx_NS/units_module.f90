!***********************************************************************
!
!     DRIVER: units_module
!     Copyright: Matthias Liebendoerfer, 30.09.2011, GNU public license.
!
!***********************************************************************

      module units_module

      type units_in_cgs
        double precision :: e,pi,sinsqthetw     !numbers
        double precision :: fm,km,MeV,s         !units
        double precision :: c,ec,ga,G,GMs_c2km,gv,kB,mb,Ms,Rs !cgs constants
        double precision :: h,GF,me,Q           !non-cgs constants
      end type units_in_cgs

      type(units_in_cgs), parameter :: units = units_in_cgs(&
           &  2.718281828459045d0,&        !e
           &  3.1415926535898d0,&          !pi
           &  0.2325d0,&                   !sinsqthetw
           &  1.d-13,&                     !fm
           &  1.d5,&                       !km
           &  1.602192d-6,&                !MeV
           &  1.d0,&                       !s
           &  2.997924562d10,&             !c
           &  4.80324214d-10,&             !ec [sqrt(erg*cm)]
           &  1.23d0,&                     !ga
           &  6.672387286039493d-08,&      !G [cm^3/s^2/g]
           &  1.47664d0,&                  !GMs_c2km
           &  1.d0,&                       !gv
           &  1.38066244d-16,&             !kB
           &  1.674d-24,&                  !mb
           &  1.989d33,&                   !Ms
           &  6.9599d10,&                  !Rs
           &  4.1356943d-21,&              !h [MeV*s]
           &  8.957d-44,&                  !GF [MeV*cm^3]
           &  0.511d0,&                    !me [MeV]
           &  1.2935d0)                    !Q=neutron-proton mass [MeV]

      end module units_module

!***********************************************************************

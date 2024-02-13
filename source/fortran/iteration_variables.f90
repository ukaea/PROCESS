module define_iteration_variables
   !! Module to define iteration variables

#ifndef dp
  use, intrinsic :: iso_fortran_env, only: dp=>real64
#endif
   implicit none

   public

   real(dp) :: DUMMY

contains

  subroutine init_define_iteration_variables
    !! Initialise module variables
    implicit none

    DUMMY = 0.0D0
  end subroutine init_define_iteration_variables

  !! lablxc(ipnvars) : labels describing iteration variables<UL>

  !---------------------------------

  subroutine init_itv_1
    !! <LI> ( 1) aspect
    use numerics, only: lablxc, boundl, boundu
    implicit none
    lablxc(1) = 'aspect        '
    boundl(1) = 1.100D0
    boundu(1) = 10.00D0
  end subroutine init_itv_1

  real(kind(1.d0)) function itv_1()
    use physics_variables, only: aspect
    implicit none
    itv_1 = aspect
  end function itv_1

  subroutine set_itv_1(ratio)
    use physics_variables, only: aspect
    implicit none
    real(kind(1.d0)) :: ratio
    aspect = ratio
  end subroutine set_itv_1

  !---------------------------------

  subroutine init_itv_2
    !! <LI> ( 2) bt
    use numerics, only: lablxc, boundl, boundu
    implicit none
    lablxc(2) = 'bt            '
    boundl(2) = 0.010D0
    boundu(2) = 30.00D0
  end subroutine init_itv_2

  real(kind(1.d0)) function itv_2()
    use physics_variables, only: bt
    implicit none
    itv_2 = bt
  end function itv_2

  subroutine set_itv_2(ratio)
    use physics_variables, only: bt
    implicit none
    real(kind(1.d0)) :: ratio
    bt = ratio
  end subroutine set_itv_2

  !---------------------------------

  subroutine init_itv_3
    !! <LI> ( 3) rmajor
    use numerics, only: lablxc, boundl, boundu
    implicit none
    lablxc(3) = 'rmajor        '
    boundl(3) = 0.100D0
    boundu(3) = 50.00D0
  end subroutine init_itv_3

  real(kind(1.d0)) function itv_3()
    use physics_variables, only: rmajor
    implicit none
    itv_3 = rmajor
  end function itv_3

  subroutine set_itv_3(ratio)
    use physics_variables, only: rmajor
    implicit none
    real(kind(1.d0)) :: ratio
    rmajor = ratio
  end subroutine set_itv_3

  !---------------------------------

  subroutine init_itv_4
    !! <LI> ( 4) te
    use numerics, only: lablxc, boundl, boundu
    implicit none
    lablxc(4) = 'te            '
    boundl(4) = 5.000D0
    boundu(4) = 150.0D0
  end subroutine init_itv_4

  real(kind(1.d0)) function itv_4()
    use physics_variables, only: te
    implicit none
    itv_4 = te
  end function itv_4

  subroutine set_itv_4(ratio)
    use physics_variables, only: te
    implicit none
    real(kind(1.d0)) :: ratio
    te = ratio
  end subroutine

  !---------------------------------

  subroutine init_itv_5
    !! <LI> ( 5) beta
    use numerics, only: lablxc, boundl, boundu
    implicit none
    lablxc(5) = 'beta          '
    boundl(5) = 0.001D0
    boundu(5) = 1.000D0
  end subroutine init_itv_5

  real(kind(1.d0)) function itv_5()
    use physics_variables, only: beta
    implicit none
    itv_5 = beta
  end function itv_5

  subroutine set_itv_5(ratio)
    use physics_variables, only: beta
    implicit none
    real(kind(1.d0)) :: ratio
    beta = ratio
  end subroutine

  !---------------------------------

  subroutine init_itv_6
    !! <LI> ( 6) dene
    use numerics, only: lablxc, boundl, boundu
    implicit none
    lablxc(6) = 'dene          '
    boundl(6) = 1.00D19
    boundu(6) = 1.00D21
  end subroutine init_itv_6

  real(kind(1.d0)) function itv_6()
    use physics_variables, only: dene
    implicit none
    itv_6 = dene
  end function itv_6

  subroutine set_itv_6(ratio)
    use physics_variables, only: dene
    implicit none
    real(kind(1.d0)) :: ratio
    dene = ratio
  end subroutine set_itv_6

  !---------------------------------

  subroutine init_itv_7
    !! <LI> ( 7) rnbeam
    use numerics, only: lablxc, boundl, boundu
    implicit none
    lablxc(7) = 'rnbeam        '
    boundl(7) = 1.00D-6
    boundu(7) = 1.000D0
  end subroutine init_itv_7

  real(kind(1.d0)) function itv_7()
    use physics_variables, only: rnbeam
    implicit none
    itv_7 = rnbeam
  end function itv_7

  subroutine set_itv_7(ratio)
    use physics_variables, only: rnbeam
    implicit none
    real(kind(1.d0)) :: ratio
    rnbeam = ratio
  end subroutine set_itv_7

  !---------------------------------

  subroutine init_itv_8
    !! <LI> ( 8) fbeta (f-value for equation 6)
    use numerics, only: lablxc, boundl, boundu
    implicit none
    lablxc(8) = 'fbeta         '
    boundl(8) = 0.001D0
    boundu(8) = 1.000D0
  end subroutine init_itv_8

  real(kind(1.d0)) function itv_8()
    use constraint_variables, only: fbeta
    implicit none
    itv_8 = fbeta
  end function itv_8

  subroutine set_itv_8(ratio)
    use constraint_variables, only: fbeta
    implicit none
    real(kind(1.d0)) :: ratio
    fbeta = ratio
  end subroutine set_itv_8

  !---------------------------------

  subroutine init_itv_9
    !! <LI> ( 9) fdene (f-value for equation 5)
    use numerics, only: lablxc, boundl, boundu
    implicit none
    lablxc(9) = 'fdene         '
    boundl(9) = 0.001D0
    boundu(9) = 1.000D0
  end subroutine init_itv_9

  real(kind(1.d0)) function itv_9()
    use constraint_variables, only: fdene
    implicit none
    itv_9 = fdene
  end function itv_9

  subroutine set_itv_9(ratio)
    use constraint_variables, only: fdene
    implicit none
    real(kind(1.d0)) :: ratio
    fdene = ratio
  end subroutine set_itv_9

  !---------------------------------

  subroutine init_itv_10
    !! <LI> (10) hfact
    use numerics, only: lablxc, boundl, boundu
    implicit none
    lablxc(10) = 'hfact         '
    boundl(10) = 0.100D0
    boundu(10) = 3.000D0
  end subroutine init_itv_10

  real(kind(1.d0)) function itv_10()
    use physics_variables, only: hfact
    implicit none
    itv_10 = hfact
  end function itv_10

  subroutine set_itv_10(ratio)
    use physics_variables, only: hfact
    implicit none
    real(kind(1.d0)) :: ratio
    hfact = ratio
  end subroutine set_itv_10

  !---------------------------------

  subroutine init_itv_11
    !! <LI> (11) pheat
    use numerics, only: lablxc, boundl, boundu
    implicit none
    lablxc(11) = 'pheat         '
    boundl(11) = 1.00D-3
    boundu(11) = 1.000D3
  end subroutine init_itv_11

  real(kind(1.d0)) function itv_11()
    use current_drive_variables, only: pheat
    implicit none
    itv_11 = pheat
  end function itv_11

  subroutine set_itv_11(ratio)
    use current_drive_variables, only: pheat
    implicit none
    real(kind(1.d0)) :: ratio
    pheat = ratio
  end subroutine set_itv_11

  !---------------------------------

  subroutine init_itv_12
    !! <LI> (12) oacdcp
    use numerics, only: lablxc, boundl, boundu
    implicit none
    lablxc(12) = 'oacdcp        '
    boundl(12) = 1.000D5
    boundu(12) =  1.500D8
  end subroutine init_itv_12

  real(kind(1.d0)) function itv_12()
    use tfcoil_variables, only: oacdcp
    implicit none
    itv_12 = oacdcp
  end function itv_12

  subroutine set_itv_12(ratio)
    use tfcoil_variables, only: oacdcp
    implicit none
    real(kind(1.d0)) :: ratio
    oacdcp = ratio
  end subroutine set_itv_12

  !---------------------------------

  subroutine init_itv_13
    !! <LI> (13) tfcth (NOT RECOMMENDED)
    use numerics, only: lablxc, boundl, boundu
    implicit none
    lablxc(13) = 'tfcth         '
    boundl(13) = 0.100D0
    boundu(13) = 5.000D0
  end subroutine init_itv_13

  real(kind(1.d0)) function itv_13()
    use build_variables, only: tfcth
    use error_handling, only: report_error
    use stellarator_variables, only: istell
    implicit none
    itv_13 = tfcth
    if (istell /= 0) then
        call report_error(46)
    end if
  end function itv_13

  subroutine set_itv_13(ratio)
    use build_variables, only: tfcth
    implicit none
    real(kind(1.d0)) :: ratio
    tfcth = ratio
  end subroutine set_itv_13

  !---------------------------------

  subroutine init_itv_14
    !! <LI> (14) fwalld (f-value for equation 8)
    use numerics, only: lablxc, boundl, boundu
    implicit none
    lablxc(14) = 'fwalld        '
    boundl(14) = 0.001D0
    boundu(14) = 1.000D0
  end subroutine init_itv_14

  real(kind(1.d0)) function itv_14()
    use constraint_variables, only: fwalld
    implicit none
    itv_14 = fwalld
  end function itv_14

  subroutine set_itv_14(ratio)
    use constraint_variables, only: fwalld
    implicit none
    real(kind(1.d0)) :: ratio
    fwalld = ratio
  end subroutine set_itv_14

  !---------------------------------

  subroutine init_itv_15
    !! <LI> (15) fvs (f-value for equation 12)
    use numerics, only: lablxc, boundl, boundu
    implicit none
    lablxc(15) = 'fvs           '
    boundl(15) = 0.001D0
    boundu(15) = 10.000
  end subroutine init_itv_15

  real(kind(1.d0)) function itv_15()
    use constraint_variables, only: fvs
    implicit none
    itv_15 = fvs
  end function itv_15

  subroutine set_itv_15(ratio)
    use constraint_variables, only: fvs
    implicit none
    real(kind(1.d0)) :: ratio
    fvs = ratio
  end subroutine set_itv_15

  !---------------------------------

  subroutine init_itv_16
    !! <LI> (16) ohcth
    use numerics, only: lablxc, boundl, boundu
    implicit none
    lablxc(16) = 'ohcth         '
    boundl(16) = 0.010D0
    boundu(16) = 10.00D0
  end subroutine init_itv_16

  real(kind(1.d0)) function itv_16()
    use build_variables, only: ohcth
    implicit none
    itv_16 = ohcth
  end function itv_16

  subroutine set_itv_16(ratio)
    use build_variables, only: ohcth
    implicit none
    real(kind(1.d0)) :: ratio
    ohcth = ratio
  end subroutine set_itv_16

  !---------------------------------

  subroutine init_itv_17
    !! <LI> (17) tdwell
    use numerics, only: lablxc, boundl, boundu
    implicit none
    lablxc(17) = 'tdwell        '
    boundl(17) = 0.100D0
    boundu(17) = 1.000D8
  end subroutine init_itv_17

  real(kind(1.d0)) function itv_17()
    use times_variables, only: tdwell
    implicit none
    itv_17 = tdwell
  end function itv_17

  subroutine set_itv_17(ratio)
    use times_variables, only: tdwell
    implicit none
    real(kind(1.d0)) :: ratio
    tdwell = ratio
  end subroutine set_itv_17

  !---------------------------------

  subroutine init_itv_18
    !! <LI> (18) q
    use numerics, only: lablxc, boundl, boundu
    implicit none
    lablxc(18) = 'q             '
    boundl(18) = 2.000D0
    boundu(18) = 50.00D0
  end subroutine init_itv_18

  real(kind(1.d0)) function itv_18()
    use physics_variables, only: q
    implicit none
    itv_18 = q
  end function itv_18

  subroutine set_itv_18(ratio)
    use physics_variables, only: q
    implicit none
    real(kind(1.d0)) :: ratio
    q = ratio
  end subroutine set_itv_18

  !---------------------------------

  subroutine init_itv_19
    !! <LI> (19) enbeam
    use numerics, only: lablxc, boundl, boundu
    implicit none
    lablxc(19) = 'enbeam        '
    boundl(19) = 1.000D0
    boundu(19) = 1.000D6
  end subroutine init_itv_19

  real(kind(1.d0)) function itv_19()
    use current_drive_variables, only: enbeam
    implicit none
    itv_19 = enbeam
  end function itv_19

  subroutine set_itv_19(ratio)
    use current_drive_variables, only: enbeam
    implicit none
    real(kind(1.d0)) :: ratio
    enbeam = ratio
  end subroutine

  !---------------------------------

  subroutine init_itv_20
    !! <LI> (20) tcpav
    use numerics, only: lablxc, boundl, boundu
    implicit none
    lablxc(20) = 'tcpav         '
    boundl(20) = 40.00D0
    boundu(20) = 3.000D2
  end subroutine init_itv_20

  real(kind(1.d0)) function itv_20()
    use tfcoil_variables, only: tcpav
    implicit none
    itv_20 = tcpav
  end function itv_20

  subroutine set_itv_20(ratio)
    use tfcoil_variables, only: tcpav
    implicit none
    real(kind(1.d0)) :: ratio
    tcpav = ratio
  end subroutine set_itv_20

  !---------------------------------

  subroutine init_itv_21
    !! <LI> (21) ftburn (f-value for equation 13)
    use numerics, only: lablxc, boundl, boundu
    implicit none
    lablxc(21) = 'ftburn        '
    boundl(21) = 0.001D0
    boundu(21) = 1.000D0
  end subroutine init_itv_21

  real(kind(1.d0)) function itv_21()
    use constraint_variables, only: ftburn
    implicit none
    itv_21 = ftburn
  end function itv_21

  subroutine set_itv_21(ratio)
    use constraint_variables, only: ftburn
    implicit none
    real(kind(1.d0)) :: ratio
    ftburn = ratio
  end subroutine set_itv_21

  !! <LI> (22) NOT USED
  !---------------------------------

  real(kind(1.d0)) function itv_22()
    implicit none
    write(*,*) 'Iteration variable 22 is no longer in use.'
  end function itv_22

  !---------------------------------

  subroutine init_itv_23
    !! <LI> (23) fcoolcp
    use numerics, only: lablxc, boundl, boundu
    implicit none
    lablxc(23) = 'fcoolcp       '
    boundl(23) = 0.100D0
    boundu(23) = 0.500D0
  end subroutine init_itv_23

  real(kind(1.d0)) function itv_23()
    use tfcoil_variables, only: fcoolcp
    implicit none
    itv_23 = fcoolcp
  end function itv_23

  subroutine set_itv_23(ratio)
    use tfcoil_variables, only: fcoolcp
    implicit none
    real(kind(1.d0)) :: ratio
    fcoolcp = ratio
  end subroutine set_itv_23

  !! <LI> (24) NOT USED
  !---------------------------------

  real(kind(1.d0)) function itv_24()
    implicit none
    write(*,*) 'Iteration variable 24 is no longer in use.'
  end function itv_24

  !---------------------------------

  subroutine init_itv_25
    !! <LI> (25) fpnetel (f-value for equation 16)
    use numerics, only: lablxc, boundl, boundu
    implicit none
    lablxc(25) = 'fpnetel       '
    boundl(25) = 0.001D0
    boundu(25) = 1.000D0
  end subroutine init_itv_25

  real(kind(1.d0)) function itv_25()
    use constraint_variables, only: fpnetel
    implicit none
    itv_25 = fpnetel
  end function itv_25

  subroutine set_itv_25(ratio)
    use constraint_variables, only: fpnetel
    implicit none
    real(kind(1.d0)) :: ratio
    fpnetel = ratio
  end subroutine set_itv_25

  !---------------------------------

  subroutine init_itv_26
    !! <LI> (26) ffuspow (f-value for equation 9)
    use numerics, only: lablxc, boundl, boundu
    implicit none
    lablxc(26) = 'ffuspow       '
    boundl(26) = 0.001D0
    boundu(26) = 1.000D0
  end subroutine init_itv_26

  real(kind(1.d0)) function itv_26()
    use constraint_variables, only: ffuspow
    implicit none
    itv_26 = ffuspow
  end function itv_26

  subroutine set_itv_26(ratio)
    use constraint_variables, only: ffuspow
    implicit none
    real(kind(1.d0)) :: ratio
    ffuspow = ratio
  end subroutine set_itv_26

  !---------------------------------

  subroutine init_itv_27
    !! <LI> (27) fhldiv (f-value for equation 18)
    use numerics, only: lablxc, boundl, boundu
    implicit none
    lablxc(27) = 'fhldiv        '
    boundl(27) = 0.001D0
    boundu(27) = 1.000D0
  end subroutine init_itv_27

  real(kind(1.d0)) function itv_27()
    use constraint_variables, only: fhldiv
    implicit none
    itv_27 =  fhldiv
  end function itv_27

  subroutine set_itv_27(ratio)
    use constraint_variables, only: fhldiv
    implicit none
    real(kind(1.d0)) :: ratio
    fhldiv = ratio
  end subroutine set_itv_27

  !---------------------------------

  subroutine init_itv_28
    !! <LI> (28) fradpwr (f-value for equation 17), total radiation fraction
    use numerics, only: lablxc, boundl, boundu
    implicit none
    lablxc(28) = 'fradpwr       '
    boundl(28) = 0.001D0
    boundu(28) = 0.990D0
  end subroutine init_itv_28

  real(kind(1.d0)) function itv_28()
    use constraint_variables, only: fradpwr
    implicit none
    itv_28 = fradpwr
  end function itv_28

  subroutine set_itv_28(ratio)
    use constraint_variables, only: fradpwr
    implicit none
    real(kind(1.d0)) :: ratio
    fradpwr = ratio
  end subroutine set_itv_28

  !---------------------------------

  subroutine init_itv_29
    !! <LI> (29) bore
    use numerics, only: lablxc, boundl, boundu
    implicit none
    lablxc(29) = 'bore          '
    boundl(29) = 0.100D0
    boundu(29) = 10.00D0
  end subroutine init_itv_29

  real(kind(1.d0)) function itv_29()
    use build_variables, only: bore
    implicit none
    itv_29 = bore
  end function itv_29

  subroutine set_itv_29(ratio)
    use build_variables, only: bore
    implicit none
    real(kind(1.d0)) :: ratio
    bore = ratio
  end subroutine set_itv_29

  !---------------------------------

  subroutine init_itv_30
    !! <LI> (30) fmva (f-value for equation 19)
    use numerics, only: lablxc, boundl, boundu
    implicit none
    lablxc(30) = 'fmva          '
    boundl(30) = 0.010D0
    boundu(30) = 1.000D0
  end subroutine init_itv_30

  real(kind(1.d0)) function itv_30()
    use constraint_variables, only: fmva
    implicit none
    itv_30 = fmva
  end function itv_30

  subroutine set_itv_30(ratio)
    use constraint_variables, only: fmva
    implicit none
    real(kind(1.d0)) :: ratio
    fmva = ratio
  end subroutine set_itv_30

  !---------------------------------

  subroutine init_itv_31
    !! <LI> (31) gapomin
    use numerics, only: lablxc, boundl, boundu
    implicit none
    lablxc(31) = 'gapomin       '
    boundl(31) = 0.001D0
    boundu(31) = 1.000D1
  end subroutine init_itv_31

  real(kind(1.d0)) function itv_31()
    use build_variables, only: gapomin
    implicit none
    itv_31 = gapomin
  end function itv_31

  subroutine set_itv_31(ratio)
    use build_variables, only: gapomin
    implicit none
    real(kind(1.d0)) :: ratio
    gapomin = ratio
  end subroutine set_itv_31

  !---------------------------------

  subroutine init_itv_32
    !! <LI> (32) frminor (f-value for equation 21)
    use numerics, only: lablxc, boundl, boundu
    implicit none
    lablxc(32) = 'frminor       '
    boundl(32) = 0.001D0
    boundu(32) = 1.000D0
  end subroutine init_itv_32

  real(kind(1.d0)) function itv_32()
    use constraint_variables, only: frminor
    implicit none
    itv_32 = frminor
  end function itv_32

  subroutine set_itv_32(ratio)
    use constraint_variables, only: frminor
    implicit none
    real(kind(1.d0)) :: ratio
    frminor = ratio
  end subroutine set_itv_32

  !---------------------------------

  subroutine init_itv_33
    !! <LI> (33) fportsz (f-value for equation 20)
    use numerics, only: lablxc, boundl, boundu
    implicit none
    lablxc(33) = 'fportsz       '
    boundl(33) = 0.001D0
    boundu(33) = 1.000D0
  end subroutine init_itv_33

  real(kind(1.d0)) function itv_33()
    use constraint_variables, only: fportsz
    implicit none
    itv_33 = fportsz
  end function itv_33

  subroutine set_itv_33(ratio)
    use constraint_variables, only: fportsz
    implicit none
    real(kind(1.d0)) :: ratio
    fportsz = ratio
  end subroutine set_itv_33

  !---------------------------------

  subroutine init_itv_34
    !! <LI> (34) fdivcol (f-value for equation 22)
    use numerics, only: lablxc, boundl, boundu
    implicit none
    lablxc(34) = 'fdivcol       '
    boundl(34) = 0.001D0
    boundu(34) = 1.000D0
  end subroutine init_itv_34

  real(kind(1.d0)) function itv_34()
    use constraint_variables, only: fdivcol
    implicit none
    itv_34 = fdivcol
  end function itv_34

  subroutine set_itv_34(ratio)
    use constraint_variables, only: fdivcol
    implicit none
    real(kind(1.d0)) :: ratio
    fdivcol = ratio
  end subroutine set_itv_34

  !---------------------------------

  subroutine init_itv_35
    !! <LI> (35) fpeakb (f-value for equation 25)
    use numerics, only: lablxc, boundl, boundu
    implicit none
    lablxc(35) = 'fpeakb        '
    boundl(35) = 0.001D0
    boundu(35) = 1.000D0
  end subroutine init_itv_35

  real(kind(1.d0)) function itv_35()
    use constraint_variables, only: fpeakb
    implicit none
    itv_35 = fpeakb
  end function itv_35

  subroutine set_itv_35(ratio)
    use constraint_variables, only: fpeakb
    implicit none
    real(kind(1.d0)) :: ratio
    fpeakb = ratio
  end subroutine set_itv_35

  !---------------------------------

  subroutine init_itv_36
    !! <LI> (36) fbetatry (f-value for equation 24)
    use numerics, only: lablxc, boundl, boundu
    implicit none
    lablxc(36) = 'fbetatry      '
    boundl(36) = 0.001D0
    boundu(36) = 1.000D0
  end subroutine init_itv_36

  real(kind(1.d0)) function itv_36()
    use constraint_variables, only: fbetatry
    implicit none
    itv_36 = fbetatry
  end function itv_36

  subroutine set_itv_36(ratio)
    use constraint_variables, only: fbetatry
    implicit none
    real(kind(1.d0)) :: ratio
    fbetatry = ratio
  end subroutine set_itv_36

  !---------------------------------

  subroutine init_itv_37
    !! <LI> (37) coheof
    use numerics, only: lablxc, boundl, boundu
    implicit none
    lablxc(37) = 'coheof        '
    boundl(37) = 1.000D5
    boundu(37) = 1.000D8
  end subroutine init_itv_37

  real(kind(1.d0)) function itv_37()
    use pfcoil_variables, only: coheof
    implicit none
    itv_37 = coheof
  end function itv_37

  subroutine set_itv_37(ratio)
    use pfcoil_variables, only: coheof
    implicit none
    real(kind(1.d0)) :: ratio
    coheof = ratio
  end subroutine set_itv_37

  !---------------------------------

  subroutine init_itv_38
    !! <LI> (38) fjohc (f-value for equation 26)
    use numerics, only: lablxc, boundl, boundu
    implicit none
    lablxc(38) = 'fjohc         '
    boundl(38) = 0.010D0
    boundu(38) = 1.000D0
  end subroutine init_itv_38

  real(kind(1.d0)) function itv_38()
    use constraint_variables, only: fjohc
    implicit none
    itv_38 = fjohc
  end function itv_38

  subroutine set_itv_38(ratio)
    use constraint_variables, only: fjohc
    implicit none
    real(kind(1.d0)) :: ratio
    fjohc = ratio
  end subroutine set_itv_38

  !---------------------------------

  subroutine init_itv_39
    !! <LI> (39) fjohc0 (f-value for equation 27)
    use numerics, only: lablxc, boundl, boundu
    implicit none
    lablxc(39) = 'fjohc0        '
    boundl(39) = 0.001D0
    boundu(39) = 1.000D0
  end subroutine init_itv_39

  real(kind(1.d0)) function itv_39()
    use constraint_variables, only: fjohc0
    implicit none
    itv_39 = fjohc0
  end function itv_39

  subroutine set_itv_39(ratio)
    use constraint_variables, only: fjohc0
    implicit none
    real(kind(1.d0)) :: ratio
    fjohc0 = ratio
  end subroutine set_itv_39

  !---------------------------------

  subroutine init_itv_40
    !! <LI> (40) fgamcd (f-value for equation 37)
    use numerics, only: lablxc, boundl, boundu
    implicit none
    lablxc(40) = 'fgamcd        '
    boundl(40) = 0.001D0
    boundu(40) = 1.000D0
  end subroutine init_itv_40

  real(kind(1.d0)) function itv_40()
    use constraint_variables, only: fgamcd
    implicit none
    itv_40 = fgamcd
  end function itv_40

  subroutine set_itv_40(ratio)
    use constraint_variables, only: fgamcd
    implicit none
    real(kind(1.d0)) :: ratio
    fgamcd = ratio
  end subroutine set_itv_40

  !---------------------------------

  subroutine init_itv_41
    !! <LI> (41) fcohbop
    use numerics, only: lablxc, boundl, boundu
    implicit none
    lablxc(41) = 'fcohbop       '
    boundl(41) = 0.001D0
    boundu(41) = 1.000D0
  end subroutine init_itv_41

  real(kind(1.d0)) function itv_41()
    use pfcoil_variables, only: fcohbop
    implicit none
    itv_41 = fcohbop
  end function itv_41

  subroutine set_itv_41(ratio)
    use pfcoil_variables, only: fcohbop
    implicit none
    real(kind(1.d0)) :: ratio
    fcohbop = ratio
  end subroutine set_itv_41

  !---------------------------------

  subroutine init_itv_42
    !! <LI> (42) gapoh
    use numerics, only: lablxc, boundl, boundu
    implicit none
    lablxc(42) = 'gapoh         '
    boundl(42) = 0.001D0
    boundu(42) = 10.00D0
  end subroutine init_itv_42

  real(kind(1.d0)) function itv_42()
    use build_variables, only: gapoh
    implicit none
    itv_42 = gapoh
  end function

  subroutine set_itv_42(ratio)
    use build_variables, only: gapoh
    implicit none
    real(kind(1.d0)) :: ratio
    gapoh = ratio
  end subroutine set_itv_42

  !! <LI> (43) NOT USED
  !---------------------------------

  real(kind(1.d0)) function itv_43()
    implicit none
    write(*,*) 'Iteration variable 43 is no longer in use.'
  end function itv_43

  !---------------------------------

  subroutine init_itv_44
    !! <LI> (44) fvsbrnni
    use numerics, only: lablxc, boundl, boundu
    implicit none
    lablxc(44) = 'fvsbrnni      '
    boundl(44) = 0.001D0
    boundu(44) = 1.000D0
  end subroutine init_itv_44

  real(kind(1.d0)) function itv_44()
    use physics_variables, only: fvsbrnni
    implicit none
    itv_44 = fvsbrnni
  end function itv_44

  subroutine set_itv_44(ratio)
    use physics_variables, only: fvsbrnni
    implicit none
    real(kind(1.d0)) :: ratio
    fvsbrnni = ratio
  end subroutine set_itv_44

  !---------------------------------

  subroutine init_itv_45
    !! <LI> (45) fqval (f-value for equation 28)
    use numerics, only: lablxc, boundl, boundu
    implicit none
    lablxc(45) = 'fqval         '
    boundl(45) = 0.001D0
    boundu(45) = 1.000D0
  end subroutine init_itv_45

  real(kind(1.d0)) function itv_45()
    use constraint_variables, only: fqval
    implicit none
    itv_45 = fqval
  end function itv_45

  subroutine set_itv_45(ratio)
    use constraint_variables, only: fqval
    implicit none
    real(kind(1.d0)) :: ratio
    fqval = ratio
  end subroutine set_itv_45

  !---------------------------------

  subroutine init_itv_46
    !! <LI> (46) fpinj (f-value for equation 30)
    use numerics, only: lablxc, boundl, boundu
    implicit none
    lablxc(46) = 'fpinj         '
    boundl(46) = 0.001D0
    boundu(46) = 1.000D0
  end subroutine init_itv_46

  real(kind(1.d0)) function itv_46()
    use constraint_variables, only: fpinj
    implicit none
    itv_46 = fpinj
  end function

  subroutine set_itv_46(ratio)
    use constraint_variables, only: fpinj
    implicit none
    real(kind(1.d0)) :: ratio
    fpinj = ratio
  end subroutine set_itv_46

  !---------------------------------

  subroutine init_itv_47
    !! <LI> (47) feffcd
    use numerics, only: lablxc, boundl, boundu
    implicit none
      lablxc(47) = 'feffcd        '
    boundl(47) = 0.001D0
    boundu(47) = 1.000D0
  end subroutine init_itv_47

  real(kind(1.d0)) function itv_47()
    use current_drive_variables, only: feffcd
    implicit none
    itv_47 = feffcd
  end function itv_47

  subroutine set_itv_47(ratio)
    use current_drive_variables, only: feffcd
    implicit none
    real(kind(1.d0)) :: ratio
    feffcd = ratio
  end subroutine set_itv_47

  !---------------------------------

  subroutine init_itv_48
    !! <LI> (48) fstrcase (f-value for equation 31)
    use numerics, only: lablxc, boundl, boundu
    implicit none
    lablxc(48) = 'fstrcase      '
    boundl(48) = 0.001D0
    boundu(48) = 1.000D0
  end subroutine init_itv_48

  real(kind(1.d0)) function itv_48()
    use constraint_variables, only: fstrcase
    implicit none
    itv_48 = fstrcase
  end function itv_48

  subroutine set_itv_48(ratio)
    use constraint_variables, only: fstrcase
    implicit none
    real(kind(1.d0)) :: ratio
    fstrcase = ratio
  end subroutine set_itv_48

  !---------------------------------

  subroutine init_itv_49
    !! <LI> (49) fstrcond (f-value for equation 32)
    use numerics, only: lablxc, boundl, boundu
    implicit none
    lablxc(49) = 'fstrcond      '
    boundl(49) = 0.001D0
    boundu(49) = 1.000D0
  end subroutine init_itv_49

  real(kind(1.d0)) function itv_49()
    use constraint_variables, only: fstrcond
    implicit none
    itv_49 = fstrcond
  end function itv_49

  subroutine set_itv_49(ratio)
    use constraint_variables, only: fstrcond
    implicit none
    real(kind(1.d0)) :: ratio
    fstrcond = ratio
  end subroutine set_itv_49

  !---------------------------------

  subroutine init_itv_50
    !! <LI> (50) fiooic (f-value for equation 33)
    use numerics, only: lablxc, boundl, boundu
    implicit none
    lablxc(50) = 'fiooic        '
    boundl(50) = 0.001D0
    boundu(50) = 1.000D0
  end subroutine init_itv_50

  real(kind(1.d0)) function itv_50()
    use constraint_variables, only: fiooic
    implicit none
    itv_50 = fiooic
  end function itv_50

  subroutine set_itv_50(ratio)
    use constraint_variables, only: fiooic
    implicit none
    real(kind(1.d0)) :: ratio
    fiooic = ratio
  end subroutine set_itv_50

  !---------------------------------

  subroutine init_itv_51
    !! <LI> (51) fvdump (f-value for equation 34)
    use numerics, only: lablxc, boundl, boundu
    implicit none
    lablxc(51) = 'fvdump        '
    boundl(51) = 0.001D0
    boundu(51) = 1.000D0
  end subroutine init_itv_51

  real(kind(1.d0)) function itv_51()
    use constraint_variables, only: fvdump
    implicit none
    itv_51 = fvdump
  end function

  subroutine set_itv_51(ratio)
    use constraint_variables, only: fvdump
    implicit none
    real(kind(1.d0)) :: ratio
    fvdump = ratio
  end subroutine set_itv_51

  !---------------------------------

  subroutine init_itv_52
    !! <LI> (52) vdalw
    use numerics, only: lablxc, boundl, boundu
    implicit none
    lablxc(52) = 'vdalw         '
    boundl(52) = 0.001D0
    boundu(52) = 1.000D6
  end subroutine init_itv_52

  real(kind(1.d0)) function itv_52()
    use tfcoil_variables, only: vdalw
    implicit none
    itv_52 =  vdalw
  end function itv_52

  subroutine set_itv_52(ratio)
    use tfcoil_variables, only: vdalw
    implicit none
    real(kind(1.d0)) :: ratio
    vdalw = ratio
  end subroutine set_itv_52

  !---------------------------------

  subroutine init_itv_53
    !! <LI> (53) fjprot (f-value for equation 35)
    use numerics, only: lablxc, boundl, boundu
    implicit none
    lablxc(53) = 'fjprot        '
    boundl(53) = 0.001D0
    boundu(53) = 1.000D0
  end subroutine init_itv_53

  real(kind(1.d0)) function itv_53()
    use constraint_variables, only: fjprot
    implicit none
    itv_53 = fjprot
  end function itv_53

  subroutine set_itv_53(ratio)
    use constraint_variables, only: fjprot
    implicit none
    real(kind(1.d0)) :: ratio
    fjprot = ratio
  end subroutine

  !---------------------------------

  subroutine init_itv_54
    !! <LI> (54) ftmargtf (f-value for equation 36)
    use numerics, only: lablxc, boundl, boundu
    implicit none
    lablxc(54) = 'ftmargtf      '
    boundl(54) = 0.001D0
    boundu(54) = 1.000D0
  end subroutine init_itv_54

  real(kind(1.d0)) function itv_54()
    use constraint_variables, only: ftmargtf
    implicit none
    itv_54 = ftmargtf
  end function itv_54

  subroutine set_itv_54(ratio)
    use constraint_variables, only: ftmargtf
    implicit none
    real(kind(1.d0)) :: ratio
    ftmargtf = ratio
  end subroutine set_itv_54

  !! <LI> (55) NOT USED
  !---------------------------------

  real(kind(1.d0)) function itv_55()
    implicit none
    write(*,*) 'Iteration variable 55 is no longer in use.'
  end function itv_55

  !---------------------------------

  subroutine init_itv_56
    !! <LI> (56) tdmptf
    use numerics, only: lablxc, boundl, boundu
    implicit none
    lablxc(56) = 'tdmptf        '
    boundl(56) = 0.100D0
    boundu(56) = 100.0D0
  end subroutine init_itv_56

  real(kind(1.d0)) function itv_56()
    use tfcoil_variables, only: tdmptf
    implicit none
    itv_56 = tdmptf
  end function itv_56

  subroutine set_itv_56(ratio)
    use tfcoil_variables, only: tdmptf
    implicit none
    real(kind(1.d0)) :: ratio
    tdmptf = ratio
  end subroutine set_itv_56

  !---------------------------------

  subroutine init_itv_57
    !! <LI> (57) thkcas
    use numerics, only: lablxc, boundl, boundu
    implicit none
    lablxc(57) = 'thkcas        '
    boundl(57) = 0.050D0
    boundu(57) = 1.000D0
  end subroutine init_itv_57

  real(kind(1.d0)) function itv_57()
    use error_handling, only: report_error
    use stellarator_variables, only: istell
    use tfcoil_variables, only: thkcas
    implicit none
    itv_57 = thkcas
    if (istell /= 0) then
        call report_error(48)
    end if
  end function itv_57

  subroutine set_itv_57(ratio)
    use tfcoil_variables, only: thkcas
    implicit none
    real(kind(1.d0)) :: ratio
    thkcas = ratio
  end subroutine

  !---------------------------------

  subroutine init_itv_58
    !! <LI> (58) thwcndut
    use numerics, only: lablxc, boundl, boundu
    implicit none
    lablxc(58) = 'thwcndut      '
    boundl(58) = 0.001D0
    boundu(58) = 0.100D0
  end subroutine init_itv_58

  real(kind(1.d0)) function itv_58()
    use tfcoil_variables, only: thwcndut
    implicit none
    itv_58 = thwcndut
  end function itv_58

  subroutine set_itv_58(ratio)
    use tfcoil_variables, only: thwcndut
    implicit none
    real(kind(1.d0)) :: ratio
    thwcndut = ratio
  end subroutine set_itv_58

  !---------------------------------

  subroutine init_itv_59
    !! <LI> (59) fcutfsu
    use numerics, only: lablxc, boundl, boundu
    implicit none
    lablxc(59) = 'fcutfsu       '
    boundl(59) = 0.001D0
    boundu(59) = 1.000D0
  end subroutine init_itv_59

  real(kind(1.d0)) function itv_59()
    use tfcoil_variables, only: fcutfsu
    implicit none
    itv_59 = fcutfsu
  end function itv_59

  subroutine set_itv_59(ratio)
    use tfcoil_variables, only: fcutfsu
    implicit none
    real(kind(1.d0)) :: ratio
    fcutfsu = ratio
  end subroutine set_itv_59

  !---------------------------------

  subroutine init_itv_60
    !! <LI> (60) cpttf
    use numerics, only: lablxc, boundl, boundu
    implicit none
    lablxc(60) = 'cpttf         '
    boundl(60) = 0.001D0
    boundu(60) = 4.000D4
  end subroutine init_itv_60

  real(kind(1.d0)) function itv_60()
    use error_handling, only: report_error
    use stellarator_variables, only: istell
    use tfcoil_variables, only: i_tf_sup, cpttf
    implicit none
    itv_60 = cpttf
    if ((istell /= 0).or.(i_tf_sup /= 1)) then
        call report_error(49)
    end if
  end function itv_60

  subroutine set_itv_60(ratio)
    use tfcoil_variables, only: cpttf
    implicit none
    real(kind(1.d0)) :: ratio
    cpttf = ratio
  end subroutine set_itv_60

  !---------------------------------

  subroutine init_itv_61
    !! <LI> (61) gapds
    use numerics, only: lablxc, boundl, boundu
    implicit none
    lablxc(61) = 'gapds         '
    boundl(61) = 0.001D0
    boundu(61) = 10.00D0
  end subroutine init_itv_61

  real(kind(1.d0)) function itv_61()
    use build_variables, only: gapds
    implicit none
    itv_61 = gapds
  end function itv_61

  subroutine set_itv_61(ratio)
    use build_variables, only: gapds
    implicit none
    real(kind(1.d0)) :: ratio
    gapds = ratio
  end subroutine set_itv_61

  !---------------------------------

  subroutine init_itv_62
    !! <LI> (62) fdtmp (f-value for equation 38)
    use numerics, only: lablxc, boundl, boundu
    implicit none
    lablxc(62) = 'fdtmp         '
    boundl(62) = 0.001D0
    boundu(62) = 1.000D0
  end subroutine init_itv_62

  real(kind(1.d0)) function itv_62()
    use constraint_variables, only: fdtmp
    implicit none
    itv_62 = fdtmp
  end function itv_62

  subroutine set_itv_62(ratio)
    use constraint_variables, only: fdtmp
    implicit none
    real(kind(1.d0)) :: ratio
    fdtmp = ratio
  end subroutine set_itv_62

  !---------------------------------

  subroutine init_itv_63
    !! <LI> (63) ftpeak (f-value for equation 39)
    use numerics, only: lablxc, boundl, boundu
    implicit none
    lablxc(63) = 'ftpeak        '
    boundl(63) = 0.001D0
    boundu(63) = 1.000D0
  end subroutine init_itv_63

  real(kind(1.d0)) function itv_63()
    use constraint_variables, only: ftpeak
    implicit none
    itv_63 = ftpeak
  end function itv_63

  subroutine set_itv_63(ratio)
    use constraint_variables, only: ftpeak
    implicit none
    real(kind(1.d0)) :: ratio
    ftpeak = ratio
  end subroutine set_itv_63

  !---------------------------------

  subroutine init_itv_64
    !! <LI> (64) fauxmn (f-value for equation 40)
    use numerics, only: lablxc, boundl, boundu
    implicit none
    lablxc(64) = 'fauxmn        '
    boundl(64) = 0.001D0
    boundu(64) = 1.000D0
  end subroutine init_itv_64

  real(kind(1.d0)) function itv_64()
    use constraint_variables, only: fauxmn
    implicit none
    itv_64 = fauxmn
  end function itv_64

  subroutine set_itv_64(ratio)
    use constraint_variables, only: fauxmn
    implicit none
    real(kind(1.d0)) :: ratio
    fauxmn = ratio
  end subroutine set_itv_64

  !---------------------------------

  subroutine init_itv_65
    !! <LI> (65) tohs
    use numerics, only: lablxc, boundl, boundu
    implicit none
    lablxc(65) = 'tohs          '
    boundl(65) = 0.100D0
    boundu(65) = 1.000D3
  end subroutine init_itv_65

  real(kind(1.d0)) function itv_65()
    use error_handling, only: report_error
    use pulse_variables, only: lpulse
    use times_variables, only: tohs
    implicit none
    itv_65 = tohs
    if (lpulse /= 1) then
        call report_error(50)
    end if
  end function itv_65

  subroutine set_itv_65(ratio)
    use times_variables, only: tohs
    implicit none
    real(kind(1.d0)) :: ratio
    tohs = ratio
  end subroutine set_itv_65

  !---------------------------------

  subroutine init_itv_66
    !! <LI> (66) ftohs (f-value for equation 41)
    use numerics, only: lablxc, boundl, boundu
    implicit none
    lablxc(66) = 'ftohs         '
    boundl(66) = 0.001D0
    boundu(66) = 1.000D0
  end subroutine init_itv_66

  real(kind(1.d0)) function itv_66()
    use constraint_variables, only: ftohs
    implicit none
    itv_66 = ftohs
  end function itv_66

  subroutine set_itv_66(ratio)
    use constraint_variables, only: ftohs
    implicit none
    real(kind(1.d0)) :: ratio
    ftohs = ratio
  end subroutine set_itv_66

  !---------------------------------

  subroutine init_itv_67
    !! <LI> (67) ftcycl (f-value for equation 42)
    use numerics, only: lablxc, boundl, boundu
    implicit none
    lablxc(67) = 'ftcycl        '
    boundl(67) = 0.001D0
    boundu(67) = 1.000D0
  end subroutine init_itv_67

  real(kind(1.d0)) function itv_67()
    use constraint_variables, only: ftcycl
    implicit none
    itv_67 = ftcycl
  end function itv_67

  subroutine set_itv_67(ratio)
    use constraint_variables, only: ftcycl
    implicit none
    real(kind(1.d0)) :: ratio
    ftcycl = ratio
  end subroutine set_itv_67

  !---------------------------------

  subroutine init_itv_68
    !! <LI> (68) fptemp (f-value for equation 44)
    use numerics, only: lablxc, boundl, boundu
    implicit none
    lablxc(68) = 'fptemp        '
    boundl(68) = 0.001D0
    boundu(68) = 1.000D0
  end subroutine init_itv_68

  real(kind(1.d0)) function itv_68()
    use constraint_variables, only: fptemp
    implicit none
    itv_68 = fptemp
  end function itv_68

  subroutine set_itv_68(ratio)
    use constraint_variables, only: fptemp
    implicit none
    real(kind(1.d0)) :: ratio
    fptemp= ratio
  end subroutine set_itv_68

  !---------------------------------

  subroutine init_itv_69
    !! <LI> (69) rcool
    use numerics, only: lablxc, boundl, boundu
    implicit none
    lablxc(69) = 'rcool         '
    boundl(69) = 0.001D0
    boundu(69) = 0.010D0
  end subroutine init_itv_69

  real(kind(1.d0)) function itv_69()
    use tfcoil_variables, only: rcool
    implicit none
    itv_69 = rcool
  end function itv_69

  subroutine set_itv_69(ratio)
    use tfcoil_variables, only: rcool
    implicit none
    real(kind(1.d0)) :: ratio
    rcool = ratio
  end subroutine set_itv_69

  !---------------------------------

  subroutine init_itv_70
    !! <LI> (70) vcool
    use numerics, only: lablxc, boundl, boundu
    implicit none
    lablxc(70) = 'vcool         '
    boundl(70) = 1.000D0
    boundu(70) = 1.000D2
  end subroutine init_itv_70

  real(kind(1.d0)) function itv_70()
    use tfcoil_variables, only: vcool
    implicit none
    itv_70 = vcool
  end function itv_70

  subroutine set_itv_70(ratio)
    use tfcoil_variables, only: vcool
    implicit none
    real(kind(1.d0)) :: ratio
    vcool = ratio
  end subroutine set_itv_70

  !---------------------------------

  subroutine init_itv_71
    !! <LI> (71) fq (f-value for equation 45)
    use numerics, only: lablxc, boundl, boundu
    implicit none
    lablxc(71) = 'fq            '
    boundl(71) = 0.001D0
    boundu(71) = 1.000D0
  end subroutine init_itv_71

  real(kind(1.d0)) function itv_71()
    use constraint_variables, only: fq
    implicit none
    itv_71 = fq
  end function itv_71

  subroutine set_itv_71(ratio)
    use constraint_variables, only: fq
    implicit none
    real(kind(1.d0)) :: ratio
    fq = ratio
  end subroutine set_itv_71

  !---------------------------------

  subroutine init_itv_72
    !! <LI> (72) fipir (f-value for equation 46)
    use numerics, only: lablxc, boundl, boundu
    implicit none
    lablxc(72) = 'fipir         '
    boundl(72) = 0.001D0
    boundu(72) = 1.000D0
  end subroutine init_itv_72

  real(kind(1.d0)) function itv_72()
    use constraint_variables, only: fipir
    implicit none
    itv_72 = fipir
  end function itv_72

  subroutine set_itv_72(ratio)
    use constraint_variables, only: fipir
    implicit none
    real(kind(1.d0)) :: ratio
    fipir = ratio
  end subroutine set_itv_72

  !---------------------------------

  subroutine init_itv_73
    !! <LI> (73) scrapli
    use numerics, only: lablxc, boundl, boundu
    implicit none
    lablxc(73) = 'scrapli       '
    boundl(73) = 0.001D0
    boundu(73) = 10.00D0
  end subroutine init_itv_73

  real(kind(1.d0)) function itv_73()
    use build_variables, only: scrapli
    implicit none
    itv_73 = scrapli
  end function itv_73

  subroutine set_itv_73(ratio)
    use build_variables, only: scrapli
    implicit none
    real(kind(1.d0)) :: ratio
    scrapli = ratio
  end subroutine set_itv_73

  !---------------------------------

  subroutine init_itv_74
    !! <LI> (74) scraplo
    use numerics, only: lablxc, boundl, boundu
    implicit none
    lablxc(74) = 'scraplo       '
    boundl(74) = 0.001D0
    boundu(74) = 10.00D0
  end subroutine init_itv_74

  real(kind(1.d0)) function itv_74()
    use build_variables, only: scraplo
    implicit none
    itv_74 = scraplo
  end function itv_74

  subroutine set_itv_74(ratio)
    use build_variables, only: scraplo
    implicit none
    real(kind(1.d0)) :: ratio
    scraplo = ratio
  end subroutine set_itv_74

  !---------------------------------

  subroutine init_itv_75
    !! <LI> (75) tfootfi
    use numerics, only: lablxc, boundl, boundu
    implicit none
    lablxc(75) = 'tfootfi       '
    boundl(75) = 0.200D0
    boundu(75) = 5.000D0
  end subroutine init_itv_75

  real(kind(1.d0)) function itv_75()
    use build_variables, only: tfootfi
    implicit none
    itv_75 = tfootfi
  end function itv_75

  subroutine set_itv_75(ratio)
    use build_variables, only: tfootfi
    implicit none
    real(kind(1.d0)) :: ratio
    tfootfi = ratio
  end subroutine set_itv_75

  !! <LI> (76) NOT USED
  !---------------------------------

  real(kind(1.d0)) function itv_76()
    implicit none
    write(*,*) 'Iteration variable 76 is no longer in use.'
  end function itv_76

  !! <LI> (77) NOT USED
  !---------------------------------

  real(kind(1.d0)) function itv_77()
    implicit none
    write(*,*) 'Iteration variable 77 is no longer in use.'
  end function itv_77

  !! <LI> (78) NOT USED
  !---------------------------------

  real(kind(1.d0)) function itv_78()
    implicit none
    write(*,*) 'Iteration variable 78 is no longer in use.'
  end function  itv_78

  !---------------------------------

  subroutine init_itv_79
    !! <LI> (79) fbetap (f-value for equation 48)
    use numerics, only: lablxc, boundl, boundu
    implicit none
    lablxc(79) = 'fbetap        '
    boundl(79) = 0.001D0
    boundu(79) = 1.000D0
  end subroutine init_itv_79

  real(kind(1.d0)) function itv_79()
    use constraint_variables, only: fbetap
    implicit none
    itv_79 = fbetap
  end function itv_79

  subroutine set_itv_79(ratio)
    use constraint_variables, only: fbetap
    implicit none
    real(kind(1.d0)) :: ratio
    fbetap = ratio
  end subroutine set_itv_79

  !! <LI> (80) NOT USED
  !---------------------------------

  real(kind(1.d0)) function itv_80()
    implicit none
    write(*,*) 'Iteration variable 80 is no longer in use.'
  end function itv_80

  !---------------------------------

  subroutine init_itv_81
    !! <LI> (81) edrive
    use numerics, only: lablxc, boundl, boundu
    implicit none
    lablxc(81) = 'edrive        '
    boundl(81) = 1.000d5
    boundu(81) = 5.000d7
  end subroutine init_itv_81

  real(kind(1.d0)) function itv_81()
    use ife_variables, only: edrive
    implicit none
    itv_81 = edrive
  end function itv_81

  subroutine set_itv_81(ratio)
    use ife_variables, only: edrive
    implicit none
    real(kind(1.d0)) :: ratio
    edrive = ratio
  end subroutine set_itv_81

  !---------------------------------

  subroutine init_itv_82
    !! <LI> (82) drveff
    use numerics, only: lablxc, boundl, boundu
    implicit none
    lablxc(82) = 'drveff        '
    boundl(82) = 0.010D0
    boundu(82) = 1.000D0
  end subroutine init_itv_82

  real(kind(1.d0)) function itv_82()
    use ife_variables, only: drveff
    implicit none
    itv_82 = drveff
  end function itv_82

  subroutine set_itv_82(ratio)
    use ife_variables, only: drveff
    implicit none
    real(kind(1.d0)) :: ratio
    drveff = ratio
  end subroutine set_itv_82

  !---------------------------------

  subroutine init_itv_83
    !! <LI> (83) tgain
    use numerics, only: lablxc, boundl, boundu
    implicit none
    lablxc(83) = 'tgain         '
    boundl(83) = 1.000D0
    boundu(83) = 500.0D0
  end subroutine init_itv_83

  real(kind(1.d0)) function itv_83()
    use ife_variables, only: tgain
    implicit none
    itv_83 = tgain
  end function itv_83

  subroutine set_itv_83(ratio)
    use ife_variables, only: tgain
    implicit none
    real(kind(1.d0)) :: ratio
    tgain = ratio
  end subroutine set_itv_83

  !---------------------------------

  subroutine init_itv_84
    !! <LI> (84) chrad
    use numerics, only: lablxc, boundl, boundu
    implicit none
    lablxc(84) = 'chrad         '
    boundl(84) = 0.100D0
    boundu(84) = 20.00D0
  end subroutine init_itv_84

  real(kind(1.d0)) function itv_84()
    use ife_variables, only: chrad
    implicit none
    itv_84 = chrad
  end function itv_84

  subroutine set_itv_84(ratio)
    use ife_variables, only: chrad
    implicit none
    real(kind(1.d0)) :: ratio
    chrad = ratio
  end subroutine set_itv_84

  !---------------------------------

  subroutine init_itv_85
    !! <LI> (85) pdrive
    use numerics, only: lablxc, boundl, boundu
    implicit none
    lablxc(85) = 'pdrive        '
    boundl(85) = 1.000D6
    boundu(85) = 200.0D6
  end subroutine init_itv_85

  real(kind(1.d0)) function itv_85()
    use ife_variables, only: pdrive
    implicit none
    itv_85 = pdrive
  end function itv_85

  subroutine set_itv_85(ratio)
    use ife_variables, only: pdrive
    implicit none
    real(kind(1.d0)) :: ratio
    pdrive = ratio
  end subroutine set_itv_85

  !---------------------------------

  subroutine init_itv_86
    !! <LI> (86) frrmax (f-value for equation 50)
    use numerics, only: lablxc, boundl, boundu
    implicit none
    lablxc(86) = 'frrmax        '
    boundl(86) = 0.001D0
    boundu(86) = 1.000D0
  end subroutine init_itv_86

  real(kind(1.d0)) function itv_86()
    use ife_variables, only: frrmax
    implicit none
    itv_86 = frrmax
  end function itv_86

  subroutine set_itv_86(ratio)
    use ife_variables, only: frrmax
    implicit none
    real(kind(1.d0)) :: ratio
    frrmax = ratio
  end subroutine set_itv_86

  !! <LI> (87) NOT USED
  !---------------------------------

  !! <LI> (88) NOT USED

  !---------------------------------

  subroutine init_itv_89
    !! <LI> (89) ftbr (f-value for equation 52)
    use numerics, only: lablxc, boundl, boundu
    implicit none
    lablxc(89) = 'ftbr          '
    boundl(89) = 0.001D0
    boundu(89) = 1.000D0
  end subroutine init_itv_89

  real(kind(1.d0)) function itv_89()
    use constraint_variables, only: ftbr
    implicit none
    itv_89 = ftbr
  end function itv_89

  subroutine set_itv_89(ratio)
    use constraint_variables, only: ftbr
    implicit none
    real(kind(1.d0)) :: ratio
    ftbr = ratio
  end subroutine set_itv_89

  !---------------------------------

  subroutine init_itv_90
    !! <LI> (90) blbuith
    use numerics, only: lablxc, boundl, boundu
    implicit none
    lablxc(90) = 'blbuith       '
    boundl(90) = 0.001D0
    boundu(90) = 2.000D0
  end subroutine init_itv_90

  real(kind(1.d0)) function itv_90()
    use build_variables, only: blbuith
    implicit none
    itv_90 = blbuith
  end function itv_90

  subroutine set_itv_90(ratio)
    use build_variables, only: blbuith
    implicit none
    real(kind(1.d0)) :: ratio
    blbuith = ratio
  end subroutine set_itv_90

  !---------------------------------

  subroutine init_itv_91
    !! <LI> (91) blbuoth
    use numerics, only: lablxc, boundl, boundu
    implicit none
    lablxc(91) = 'blbuoth       '
    boundl(91) = 0.001D0
    boundu(91) = 2.000D0
  end subroutine init_itv_91

  real(kind(1.d0)) function itv_91()
    use build_variables, only: blbuoth
    implicit none
    itv_91 = blbuoth
  end function itv_91

  subroutine set_itv_91(ratio)
    use build_variables, only: blbuoth
    implicit none
    real(kind(1.d0)) :: ratio
    blbuoth = ratio
  end subroutine set_itv_91

  !---------------------------------

  subroutine init_itv_92
    !! <LI> (92) fflutf (f-value for equation 53)
    use numerics, only: lablxc, boundl, boundu
    implicit none
    lablxc(92) = 'fflutf        '
    boundl(92) = 0.001D0
    boundu(92) = 1.000D0
  end subroutine init_itv_92

  real(kind(1.d0)) function itv_92()
    use constraint_variables, only: fflutf
    implicit none
    itv_92 = fflutf
  end function itv_92

  subroutine set_itv_92(ratio)
    use constraint_variables, only: fflutf
    implicit none
    real(kind(1.d0)) :: ratio
    fflutf = ratio
  end subroutine set_itv_92

  !---------------------------------

  subroutine init_itv_93
    !! <LI> (93) shldith
    use numerics, only: lablxc, boundl, boundu
    implicit none
    lablxc(93) = 'shldith       '
    boundl(93) = 0.001D0
    boundu(93) = 10.00D0
  end subroutine init_itv_93

  real(kind(1.d0)) function itv_93()
    use build_variables, only: shldith
    implicit none
    itv_93 = shldith
  end function itv_93

  subroutine set_itv_93(ratio)
    use build_variables, only: shldith
    implicit none
    real(kind(1.d0)) :: ratio
    shldith = ratio
  end subroutine set_itv_93

  !---------------------------------

  subroutine init_itv_94
    !! <LI> (94) shldoth
    use numerics, only: lablxc, boundl, boundu
    implicit none
    lablxc(94) = 'shldoth       '
    boundl(94) = 0.001D0
    boundu(94) = 10.00D0
  end subroutine init_itv_94

  real(kind(1.d0)) function itv_94()
    use build_variables, only: shldoth
    implicit none
    itv_94 = shldoth
  end function itv_94

  subroutine set_itv_94(ratio)
    use build_variables, only: shldoth
    implicit none
    real(kind(1.d0)) :: ratio
    shldoth = ratio
  end subroutine set_itv_94

  !---------------------------------

  subroutine init_itv_95
    !! <LI> (95) fptfnuc (f-value for equation 54)
    use numerics, only: lablxc, boundl, boundu
    implicit none
    lablxc(95) = 'fptfnuc       '
    boundl(95) = 0.001D0
    boundu(95) = 1.000D0
  end subroutine init_itv_95

  real(kind(1.d0)) function itv_95()
    use constraint_variables, only: fptfnuc
    implicit none
    itv_95 = fptfnuc
  end function itv_95

  subroutine set_itv_95(ratio)
    use constraint_variables, only: fptfnuc
    implicit none
    real(kind(1.d0)) :: ratio
    fptfnuc = ratio
  end subroutine set_itv_95

  !---------------------------------

  subroutine init_itv_96
    !! <LI> (96) fvvhe (f-value for equation 55)
    use numerics, only: lablxc, boundl, boundu
    implicit none
    lablxc(96) = 'fvvhe         '
    boundl(96) = 0.001D0
    boundu(96) = 1.000D0
  end subroutine init_itv_96

  real(kind(1.d0)) function itv_96()
    use constraint_variables, only: fvvhe
    implicit none
    itv_96 = fvvhe
  end function itv_96

  subroutine set_itv_96(ratio)
    use constraint_variables, only: fvvhe
    implicit none
    real(kind(1.d0)) :: ratio
    fvvhe = ratio
  end subroutine set_itv_96

  !---------------------------------

  subroutine init_itv_97
    !! <LI> (97) fpsepr (f-value for equation 56)
    use numerics, only: lablxc, boundl, boundu
    implicit none
    lablxc(97) = 'fpsepr        '
    boundl(97) = 0.001D0
    boundu(97) = 1.000D0
  end subroutine init_itv_97

  real(kind(1.d0)) function itv_97()
    use constraint_variables, only: fpsepr
    implicit none
    itv_97 = fpsepr
  end function itv_97

  subroutine set_itv_97(ratio)
    use constraint_variables, only: fpsepr
    implicit none
    real(kind(1.d0)) :: ratio
    fpsepr = ratio
  end subroutine set_itv_97

  !---------------------------------

  subroutine init_itv_98
    !! <LI> (98) li6enrich
    use numerics, only: lablxc, boundl, boundu
    implicit none
      lablxc(98) = 'li6enrich     '
    boundl(98) = 10.00D0
    boundu(98) = 100.0D0
  end subroutine init_itv_98

  real(kind(1.d0)) function itv_98()
    use fwbs_variables, only: li6enrich
    implicit none
    itv_98 = li6enrich
  end function itv_98

  subroutine set_itv_98(ratio)
    use fwbs_variables, only: li6enrich
    implicit none
    real(kind(1.d0)) :: ratio
    li6enrich = ratio
  end subroutine set_itv_98

  !! <LI> (99) NOT USED
  !---------------------------------

  !! <LI> (100) NOT USED
  !! <LI> (101) NOT USED

  !---------------------------------



  !---------------------------------

  subroutine init_itv_103
    !! <LI> (103) flhthresh (f-value for equation 15)
    use numerics, only: lablxc, boundl, boundu
    implicit none
    lablxc(103) = 'flhthresh     '
    boundl(103) = 1.000D0
    boundu(103) = 1.000D6
  end subroutine init_itv_103

  real(kind(1.d0)) function itv_103()
    use constraint_variables, only: flhthresh
    implicit none
    itv_103 = flhthresh
  end function itv_103

  subroutine set_itv_103(ratio)
    use constraint_variables, only: flhthresh
    implicit none
    real(kind(1.d0)) :: ratio
    flhthresh = ratio
  end subroutine set_itv_103

  !---------------------------------

  subroutine init_itv_104
    !! <LI> (104) fcwr (f-value for equation 23)
    use numerics, only: lablxc, boundl, boundu
    implicit none
    lablxc(104) = 'fcwr          '
    boundl(104) = 0.001D0
    boundu(104) = 1.000D0
  end subroutine init_itv_104

  real(kind(1.d0)) function itv_104()
    use constraint_variables, only: fcwr
    implicit none
    itv_104 = fcwr
  end function itv_104

  subroutine set_itv_104(ratio)
    use constraint_variables, only: fcwr
    implicit none
    real(kind(1.d0)) :: ratio
    fcwr = ratio
  end subroutine set_itv_104

  !---------------------------------

  subroutine init_itv_105
    !! <LI> (105) fnbshinef (f-value for equation 59)
    use numerics, only: lablxc, boundl, boundu
    implicit none
    lablxc(105) = 'fnbshinef     '
    boundl(105) = 0.001D0
    boundu(105) = 1.000D0
  end subroutine init_itv_105

  real(kind(1.d0)) function itv_105()
    use constraint_variables, only: fnbshinef
    implicit none
    itv_105 = fnbshinef
  end function itv_105

  subroutine set_itv_105(ratio)
    use constraint_variables, only: fnbshinef
    implicit none
    real(kind(1.d0)) :: ratio
    fnbshinef = ratio
  end subroutine set_itv_105

  !---------------------------------

  subroutine init_itv_106
    !! <LI> (106) ftmargoh (f-value for equation 60)
    use numerics, only: lablxc, boundl, boundu
    implicit none
    lablxc(106) = 'ftmargoh      '
    boundl(106) = 0.001D0
    boundu(106) = 1.000D0
  end subroutine init_itv_106

  real(kind(1.d0)) function itv_106()
    use constraint_variables, only: ftmargoh
    implicit none
    itv_106 = ftmargoh
  end function itv_106

  subroutine set_itv_106(ratio)
    use constraint_variables, only: ftmargoh
    implicit none
    real(kind(1.d0)) :: ratio
    ftmargoh = ratio
  end subroutine set_itv_106

  !---------------------------------

  subroutine init_itv_107
    !! <LI> (107) favail (f-value for equation 61)
    use numerics, only: lablxc, boundl, boundu
    implicit none
    lablxc(107) = 'favail        '
    boundl(107) = 0.001D0
    boundu(107) = 1.000D0
  end subroutine init_itv_107

  real(kind(1.d0)) function itv_107()
    use cost_variables, only: favail
    implicit none
    itv_107 = favail
  end function itv_107

  subroutine set_itv_107(ratio)
    use cost_variables, only: favail
    implicit none
    real(kind(1.d0)) :: ratio
    favail = ratio
  end subroutine set_itv_107

  !---------------------------------

  subroutine init_itv_108
    !! <LI> (108) breeder_f: Volume of Li4SiO4 / (Volume of Be12Ti + Li4SiO4)
    use numerics, only: lablxc, boundl, boundu
    implicit none
    lablxc(108) = 'breeder_f     '
    boundl(108) = 0.060D0
    boundu(108) = 1.000D0
  end subroutine init_itv_108

  real(kind(1.d0)) function itv_108()
    use fwbs_variables, only: breeder_f
    implicit none
    itv_108 = breeder_f
  end function itv_108

  subroutine set_itv_108(ratio)
    use fwbs_variables, only: breeder_f
    implicit none
    real(kind(1.d0)) :: ratio
    breeder_f = ratio
  end subroutine set_itv_108

  !---------------------------------

  subroutine init_itv_109
    !! <LI> (109) ralpne: thermal alpha density / electron density
    use numerics, only: lablxc, boundl, boundu
    implicit none
    lablxc(109) = 'ralpne        '
    boundl(109) = 0.050D0
    boundu(109) = 0.150D0
  end subroutine init_itv_109

  real(kind(1.d0)) function itv_109()
    use physics_variables, only: ralpne
    implicit none
    itv_109 = ralpne
  end function itv_109

  subroutine set_itv_109(ratio)
    use physics_variables, only: ralpne
    implicit none
    real(kind(1.d0)) :: ratio
    ralpne = ratio
  end subroutine set_itv_109

  !---------------------------------

  subroutine init_itv_110
    !! <LI> (110) ftaulimit: Lower limit on taup/taueff the ratio of alpha
    use numerics, only: lablxc, boundl, boundu
    implicit none
    !!      particle to energy confinement times (f-value for equation 62)
    lablxc(110) = 'ftaulimit     '
    boundl(110) = 0.001D0
    boundu(110) = 1.000D0
  end subroutine init_itv_110

  real(kind(1.d0)) function itv_110()
    use constraint_variables, only: ftaulimit
    implicit none
    itv_110 = ftaulimit
  end function itv_110

  subroutine set_itv_110(ratio)
    use constraint_variables, only: ftaulimit
    implicit none
    real(kind(1.d0)) :: ratio
    ftaulimit = ratio
  end subroutine set_itv_110

  !---------------------------------

  subroutine init_itv_111
    !! <LI> (111) fniterpump: f-value for constraint that number
    use numerics, only: lablxc, boundl, boundu
    implicit none
    !!       of vacuum pumps <  TF coils (f-value for equation 63)
    lablxc(111) = 'fniterpump    '
    boundl(111) = 0.001D0
    boundu(111) = 1.000D0
  end subroutine init_itv_111

  real(kind(1.d0)) function itv_111()
    use constraint_variables, only: fniterpump
    implicit none
    itv_111 = fniterpump
  end function itv_111

  subroutine set_itv_111(ratio)
    use constraint_variables, only: fniterpump
    implicit none
    real(kind(1.d0)) :: ratio
    fniterpump = ratio
  end subroutine set_itv_111

  !---------------------------------

  subroutine init_itv_112
    !! <LI> (112) fzeffmax: f-value for max Zeff (f-value for equation 64)
    use numerics, only: lablxc, boundl, boundu
    implicit none
    lablxc(112) = 'fzeffmax      '
    boundl(112) = 0.001D0
    boundu(112) = 1.000D0
  end subroutine init_itv_112

  real(kind(1.d0)) function itv_112()
    use constraint_variables, only: fzeffmax
    implicit none
    itv_112 =  fzeffmax
  end function itv_112

  subroutine set_itv_112(ratio)
    use constraint_variables, only: fzeffmax
    implicit none
    real(kind(1.d0)) :: ratio
    fzeffmax = ratio
  end subroutine set_itv_112

  !---------------------------------

  subroutine init_itv_113
    !! <LI> (113) ftaucq: f-value for maximum VV stress (f-value for equation 65)
    use numerics, only: lablxc, boundl, boundu
    implicit none
    lablxc(113) = 'fmaxvvstress  '
    boundl(113) = 0.001D0
    boundu(113) = 1.000D0
  end subroutine init_itv_113

  real(kind(1.d0)) function itv_113()
    use constraint_variables, only: fmaxvvstress
    implicit none
    itv_113 = fmaxvvstress
  end function itv_113

  subroutine set_itv_113(ratio)
    use constraint_variables, only: fmaxvvstress
    implicit none
    real(kind(1.d0)) :: ratio
    fmaxvvstress = ratio
  end subroutine set_itv_113

  !---------------------------------

  subroutine init_itv_114
    !! <LI> (114) fw_channel_length: Length of a single first wall channel
    use numerics, only: lablxc, boundl, boundu
    implicit none
    lablxc(114) = 'fw_channel_length  '
    boundl(114) = 0.001D0
    boundu(114) = 1.000D3
  end subroutine init_itv_114

  real(kind(1.d0)) function itv_114()
    use fwbs_variables, only: fw_channel_length
    implicit none
    itv_114 = fw_channel_length
  end function itv_114

  subroutine set_itv_114(ratio)
    use fwbs_variables, only: fw_channel_length
    implicit none
    real(kind(1.d0)) :: ratio
    fw_channel_length = ratio
  end subroutine set_itv_114

  !---------------------------------

  subroutine init_itv_115
    !! <LI> (115) fpoloidalpower: f-value for max rate of change of
    use numerics, only: lablxc, boundl, boundu
    implicit none
    !!            energy in poloidal field (f-value for equation 66)
    lablxc(115) = 'fpoloidalpower'
    boundl(115) = 0.001D0
    boundu(115) = 1.000D0
  end subroutine init_itv_115

  real(kind(1.d0)) function itv_115()
    use constraint_variables, only: fpoloidalpower
    implicit none
    itv_115 = fpoloidalpower
  end function itv_115

  subroutine set_itv_115(ratio)
    use constraint_variables, only: fpoloidalpower
    implicit none
    real(kind(1.d0)) :: ratio
    fpoloidalpower = ratio
  end subroutine set_itv_115

  !---------------------------------

  subroutine init_itv_116
    !! <LI> (116) fradwall: f-value for radiation wall load limit (eq. 67)
    use numerics, only: lablxc, boundl, boundu
    implicit none
    lablxc(116) = 'fradwall      '
    boundl(116) = 0.001D0
    boundu(116) = 1.000D0
  end subroutine init_itv_116

  real(kind(1.d0)) function itv_116()
    use constraint_variables, only: fradwall
    implicit none
    itv_116 = fradwall
  end function itv_116

  subroutine set_itv_116(ratio)
    use constraint_variables, only: fradwall
    implicit none
    real(kind(1.d0)) :: ratio
    fradwall = ratio
  end subroutine set_itv_116

  !---------------------------------

  subroutine init_itv_117
    !! <LI> (117) fpsepbqar: f-value for  Psep*Bt/qar upper limit (eq. 68)
    use numerics, only: lablxc, boundl, boundu
    implicit none
    lablxc(117) = 'fpsepbqar     '
    boundl(117) = 0.001D0
    boundu(117) = 1.000D0
  end subroutine init_itv_117

  real(kind(1.d0)) function itv_117()
    use constraint_variables, only: fpsepbqar
    implicit none
    itv_117 = fpsepbqar
  end function itv_117

  subroutine set_itv_117(ratio)
    use constraint_variables, only: fpsepbqar
    implicit none
    real(kind(1.d0)) :: ratio
    fpsepbqar = ratio
  end subroutine set_itv_117

  !---------------------------------

  subroutine init_itv_118
    !! <LI> (118) fpsep: f-value to ensure separatrix power is less than
    use numerics, only: lablxc, boundl, boundu
    implicit none
    !!            value from Kallenbach divertor (f-value for equation 69)
      lablxc(118) = 'fpsep         '
    boundl(118) = 0.001D0
    boundu(118) = 1.000D0
  end subroutine init_itv_118

  real(kind(1.d0)) function itv_118()
    use constraint_variables, only: fpsep
    implicit none
    itv_118 = fpsep
  end function itv_118

  subroutine set_itv_118(ratio)
    use constraint_variables, only: fpsep
    implicit none
    real(kind(1.d0)) :: ratio
    fpsep = ratio
  end subroutine set_itv_118

  !---------------------------------

  subroutine init_itv_119
    !! <LI> (119) tesep:  separatrix temperature calculated by the Kallenbach divertor model
    use numerics, only: lablxc, boundl, boundu
    implicit none
    lablxc(119) = 'tesep         '
    boundl(119) = 0.000D0
    boundu(119) = 1.000D1
  end subroutine init_itv_119

  real(kind(1.d0)) function itv_119()
    use physics_variables, only: tesep
    implicit none
    itv_119 = tesep
  end function itv_119

  subroutine set_itv_119(ratio)
    use physics_variables, only: tesep
    implicit none
  real(kind(1.d0)) :: ratio
  tesep = ratio
  end subroutine set_itv_119

  !---------------------------------

  subroutine init_itv_120
    !! <LI> (120) ttarget: Plasma temperature adjacent to divertor sheath [eV]
    use numerics, only: lablxc, boundl, boundu
    implicit none
    lablxc(120) = 'ttarget       '
    boundl(120) = 1.000D0
    boundu(120) = 1.000D4
  end subroutine init_itv_120

  real(kind(1.d0)) function itv_120()
    ! From Kallenbach model, should be reserved if the model is going to be added back
  end function itv_120

  subroutine set_itv_120(ratio)
    implicit none
    real(kind(1.d0)) :: ratio
    ! From Kallenbach model, should be reserved if the model is going to be added back
  end subroutine set_itv_120

  !---------------------------------

  subroutine init_itv_121
    !! <LI> (121) neratio: ratio of mean SOL density at OMP to separatrix density at OMP
    use numerics, only: lablxc, boundl, boundu
    implicit none
    lablxc(121) = 'neratio       '
    boundl(121) = 0.001D0
    boundu(121) = 1.000D0
  end subroutine init_itv_121

  real(kind(1.d0)) function itv_121()
    ! From Kallenbach model, should be reserved if the model is going to be added back
  end function itv_121

  subroutine set_itv_121(ratio)
    implicit none
    real(kind(1.d0)) :: ratio
    ! From Kallenbach model, should be reserved if the model is going to be added back
  end subroutine set_itv_121

  !---------------------------------

  subroutine init_itv_122
    !! <LI> (122) oh_steel_frac : streel fraction of Central Solenoid
    use numerics, only: lablxc, boundl, boundu
    implicit none
    lablxc(122) = 'oh_steel_frac '
    boundl(122) = 0.001D0
    boundu(122) = 0.950D0
  end subroutine init_itv_122

  real(kind(1.d0)) function itv_122()
    use pfcoil_variables, only: oh_steel_frac
    implicit none
    itv_122 = oh_steel_frac
  end function itv_122

  subroutine set_itv_122(ratio)
    use pfcoil_variables, only: oh_steel_frac
    implicit none
    real(kind(1.d0)) :: ratio
    oh_steel_frac = ratio
  end subroutine set_itv_122

  !---------------------------------

  subroutine init_itv_123
    !! <LI> (123) foh_stress : f-value for CS coil Tresca yield criterion (f-value for eq. 72)
    use numerics, only: lablxc, boundl, boundu
    implicit none
    lablxc(123) = 'foh_stress    '
    boundl(123) = 0.001D0
    boundu(123) = 1.000D0
  end subroutine init_itv_123

  real(kind(1.d0)) function itv_123()
    use constraint_variables, only: foh_stress
    implicit none
    itv_123 = foh_stress
  end function itv_123

  subroutine set_itv_123(ratio)
    use constraint_variables, only: foh_stress
    implicit none
    real(kind(1.d0)) :: ratio
    foh_stress = ratio
  end subroutine set_itv_123

  !---------------------------------

  subroutine init_itv_124
    !! <LI> (124) qtargettotal : Power density on target including surface recombination [W/m2]
    use numerics, only: lablxc, boundl, boundu
    implicit none
    lablxc(124) = 'qtargettotal  '
    boundl(124) = 0.001D0
    boundu(124) = 1.000D7
  end subroutine init_itv_124

  real(kind(1.d0)) function itv_124()
    ! From Kallenbach model, should be reserved if the model is going to be added back
  end function itv_124

  subroutine set_itv_124(ratio)
    implicit none
    real(kind(1.d0)) :: ratio
    ! From Kallenbach model, should be reserved if the model is going to be added back
  end subroutine set_itv_124

  !---------------------------------

  subroutine init_itv_125
    !! <LI> (125) fimp(3) :  Beryllium density fraction relative to electron density
    use numerics, only: lablxc, boundl, boundu
    implicit none
    lablxc(125) = 'fimp(03)      '
    boundl(125) = 1.00D-8
    boundu(125) = 0.010D0
  end subroutine init_itv_125

  real(kind(1.d0)) function itv_125()
    use impurity_radiation_module, only: impurity_arr_frac
    implicit none
    itv_125 = impurity_arr_frac(3)
  end function itv_125

  subroutine set_itv_125(ratio)
    use impurity_radiation_module, only: impurity_arr_frac
    implicit none
    real(kind(1.d0)) :: ratio
    impurity_arr_frac(3) = ratio
  end subroutine set_itv_125

  !---------------------------------

  subroutine init_itv_126
    !! <LI> (126) fimp(4) :  Carbon density fraction relative to electron density
    use numerics, only: lablxc, boundl, boundu
    implicit none
    lablxc(126) = 'fimp(04)      '
    boundl(126) = 1.00D-8
    boundu(126) = 0.010D0
  end subroutine init_itv_126

  real(kind(1.d0)) function itv_126()
    use impurity_radiation_module, only: impurity_arr_frac
    implicit none
    itv_126 = impurity_arr_frac(4)
  end function itv_126

  subroutine set_itv_126(ratio)
    use impurity_radiation_module, only: impurity_arr_frac
    implicit none
    real(kind(1.d0)) :: ratio
    impurity_arr_frac(4) = ratio
  end subroutine set_itv_126

  !---------------------------------

  subroutine init_itv_127
    !! <LI> (127) fimp(5) :  Nitrogen fraction relative to electron density
    use numerics, only: lablxc, boundl, boundu
    implicit none
    lablxc(127) = 'fimp(05)      '
    boundl(127) = 1.00D-8
    boundu(127) = 0.010D0
  end subroutine init_itv_127

  real(kind(1.d0)) function itv_127()
    use impurity_radiation_module, only: impurity_arr_frac
    implicit none
    itv_127 = impurity_arr_frac(5)
  end function itv_127

  subroutine set_itv_127(ratio)
    use impurity_radiation_module, only: impurity_arr_frac
    implicit none
    real(kind(1.d0)) :: ratio
    impurity_arr_frac(5) = ratio
  end subroutine set_itv_127

  !---------------------------------

  subroutine init_itv_128
    !! <LI> (128) fimp(6) :  Oxygen density fraction relative to electron density
    use numerics, only: lablxc, boundl, boundu
    implicit none
    lablxc(128) = 'fimp(06)      '
    boundl(128) = 1.00D-8
    boundu(128) = 0.010D0
  end subroutine init_itv_128

  real(kind(1.d0)) function itv_128()
    use impurity_radiation_module, only: impurity_arr_frac
    implicit none
    itv_128 = impurity_arr_frac(6)
  end function itv_128

  subroutine set_itv_128(ratio)
    use impurity_radiation_module, only: impurity_arr_frac
    implicit none
    real(kind(1.d0)) :: ratio
    impurity_arr_frac(6) = ratio
  end subroutine set_itv_128

  !---------------------------------

  subroutine init_itv_129
    !! <LI> (129) fimp(7) :  Neon density fraction relative to electron density
    use numerics, only: lablxc, boundl, boundu
    implicit none
    lablxc(129) = 'fimp(07)      '
    boundl(129) = 1.00D-8
    boundu(129) = 0.010D0
  end subroutine init_itv_129

  real(kind(1.d0)) function itv_129()
    use impurity_radiation_module, only: impurity_arr_frac
    implicit none
    itv_129 = impurity_arr_frac(7)
  end function itv_129

  subroutine set_itv_129(ratio)
    use impurity_radiation_module, only: impurity_arr_frac
    implicit none
    real(kind(1.d0)) :: ratio
    impurity_arr_frac(7) = ratio
  end subroutine set_itv_129

  !---------------------------------

  subroutine init_itv_130
    !! <LI> (130) fimp(8) :  Silicon density fraction relative to electron density
    use numerics, only: lablxc, boundl, boundu
    implicit none
    lablxc(130) = 'fimp(08)      '
    boundl(130) = 1.00D-8
    boundu(130) = 0.010D0
  end subroutine init_itv_130

  real(kind(1.d0)) function itv_130()
    use impurity_radiation_module, only: impurity_arr_frac
    implicit none
    itv_130 = impurity_arr_frac(8)
  end function itv_130

  subroutine set_itv_130(ratio)
    use impurity_radiation_module, only: impurity_arr_frac
    implicit none
    real(kind(1.d0)) :: ratio
    impurity_arr_frac(8) = ratio
  end subroutine set_itv_130

  !---------------------------------

  subroutine init_itv_131
    !! <LI> (131) fimp(9) :  Argon density fraction relative to electron density
    use numerics, only: lablxc, boundl, boundu
    implicit none
    lablxc(131) = 'fimp(09)      '
    boundl(131) = 1.00D-8
    boundu(131) = 0.010D0
  end subroutine init_itv_131

  real(kind(1.d0)) function itv_131()
    use impurity_radiation_module, only: impurity_arr_frac
    implicit none
    itv_131 = impurity_arr_frac(9)
  end function itv_131

  subroutine set_itv_131(ratio)
    use impurity_radiation_module, only: impurity_arr_frac
    implicit none
    real(kind(1.d0)) :: ratio
    impurity_arr_frac(9) = ratio
  end subroutine set_itv_131

  !---------------------------------

  subroutine init_itv_132
    !! <LI> (132) fimp(10) :  Iron density fraction relative to electron density
    use numerics, only: lablxc, boundl, boundu
    implicit none
    lablxc(132) = 'fimp(10)      '
    boundl(132) = 1.00D-8
    boundu(132) = 0.010D0
  end subroutine init_itv_132

  real(kind(1.d0)) function itv_132()
    use impurity_radiation_module, only: impurity_arr_frac
    implicit none
    itv_132 = impurity_arr_frac(10)
  end function itv_132

  subroutine set_itv_132(ratio)
    use impurity_radiation_module, only: impurity_arr_frac
    implicit none
    real(kind(1.d0)) :: ratio
    impurity_arr_frac(10) = ratio
  end subroutine set_itv_132

  !---------------------------------

  subroutine init_itv_133
    !! <LI> (133) fimp(11) :  Nickel density fraction relative to electron density
    use numerics, only: lablxc, boundl, boundu
    implicit none
    lablxc(133) = 'fimp(11)      '
    boundl(133) = 1.00D-8
    boundu(133) = 0.010D0
  end subroutine init_itv_133

  real(kind(1.d0)) function itv_133()
    use impurity_radiation_module, only: impurity_arr_frac
    implicit none
    itv_133 = impurity_arr_frac(11)
  end function itv_133

  subroutine set_itv_133(ratio)
    use impurity_radiation_module, only: impurity_arr_frac
    implicit none
    real(kind(1.d0)) :: ratio
    impurity_arr_frac(11) = ratio
  end subroutine set_itv_133

  !---------------------------------

  subroutine init_itv_134
    !! <LI> (134) fimp(12) :  Krypton density fraction relative to electron density
    use numerics, only: lablxc, boundl, boundu
    implicit none
    lablxc(134) = 'fimp(12)      '
    boundl(134) = 1.00D-8
    boundu(134) = 0.010D0
  end subroutine init_itv_134

  real(kind(1.d0)) function itv_134()
    use impurity_radiation_module, only: impurity_arr_frac
    implicit none
    itv_134 = impurity_arr_frac(12)
  end function itv_134

  subroutine set_itv_134(ratio)
    use impurity_radiation_module, only: impurity_arr_frac
    implicit none
    real(kind(1.d0)) :: ratio
    impurity_arr_frac(12) = ratio
  end subroutine set_itv_134

  !---------------------------------

  subroutine init_itv_135
    !! <LI> (135) fimp(13) :  Xenon density fraction relative to electron density
    use numerics, only: lablxc, boundl, boundu
    implicit none
    lablxc(135) = 'fimp(13)      '
    boundl(135) = 1.00D-8
    boundu(135) = 0.010D0
  end subroutine init_itv_135

  real(kind(1.d0)) function itv_135()
    use impurity_radiation_module, only: impurity_arr_frac
    implicit none
    itv_135 = impurity_arr_frac(13)
  end function itv_135

  subroutine set_itv_135(ratio)
    use impurity_radiation_module, only: impurity_arr_frac
    implicit none
    real(kind(1.d0)) :: ratio
    impurity_arr_frac(13) = ratio
  end subroutine set_itv_135

  !---------------------------------

  subroutine init_itv_136
    !! <LI> (136) fimp(14) :  Tungsten density fraction relative to electron density
    use numerics, only: lablxc, boundl, boundu
    implicit none
    lablxc(136) = 'fimp(14)      '
    boundl(136) = 1.00D-8
    boundu(136) = 0.010D0
  end subroutine init_itv_136

  real(kind(1.d0)) function itv_136()
    use impurity_radiation_module, only: impurity_arr_frac
    implicit none
    itv_136 = impurity_arr_frac(14)
  end function itv_136

  subroutine set_itv_136(ratio)
    use impurity_radiation_module, only: impurity_arr_frac
    implicit none
    real(kind(1.d0)) :: ratio
    impurity_arr_frac(14) = ratio
  end subroutine set_itv_136

  !---------------------------------

  subroutine init_itv_137
    !! <LI> (137) fplhsep (f-value for equation 73)
    use numerics, only: lablxc, boundl, boundu
    implicit none
    lablxc(137) = 'fplhsep       '
    boundl(137) = 0.001D0
    boundu(137) = 1.000D0
  end subroutine init_itv_137

  real(kind(1.d0)) function itv_137()
    use physics_variables, only: fplhsep
    implicit none
    itv_137 = fplhsep
  end function itv_137

  subroutine set_itv_137(ratio)
    use physics_variables, only: fplhsep
    implicit none
    real(kind(1.d0)) :: ratio
    fplhsep = ratio
  end subroutine set_itv_137

  !---------------------------------

  subroutine init_itv_138
    !! <LI> (138) rebco_thickness : thickness of REBCO layer in tape (m)
    use numerics, only: lablxc, boundl, boundu
    implicit none
    lablxc(138) = 'rebco_thickness'
    boundl(138) = 0.01D-6
    boundu(138) = 100.0D-6
  end subroutine init_itv_138

  real(kind(1.d0)) function itv_138()
    use physics_variables, only: fplhsep
    implicit none
    itv_138 = fplhsep
  end function itv_138

  subroutine set_itv_138(ratio)
    use physics_variables, only: fplhsep
    implicit none
    real(kind(1.d0)) :: ratio
    fplhsep = ratio
  end subroutine set_itv_138

  !---------------------------------

  subroutine init_itv_139
    !! <LI> (139) copper_thick : thickness of copper layer in tape (m)
    use numerics, only: lablxc, boundl, boundu
    implicit none
    lablxc(139) = 'copper_thick  '
    boundl(139) = 1.00D-6
    boundu(139) = 1.00D-3
  end subroutine init_itv_139

  real(kind(1.d0)) function itv_139()
    use rebco_variables, only: copper_thick
    implicit none
    itv_139 = copper_thick
  end function itv_139

  subroutine set_itv_139(ratio)
    use rebco_variables, only: copper_thick
    implicit none
    real(kind(1.d0)) :: ratio
    copper_thick = ratio
  end subroutine set_itv_139

  !---------------------------------

  subroutine init_itv_140
    !! <LI> (140) dr_tf_wp : radial thickness of TFC winding pack (m)
    use numerics, only: lablxc, boundl, boundu
    implicit none
    lablxc(140) = 'dr_tf_wp         '
    boundl(140) = 0.001D0
    boundu(140) = 2.000D0
  end subroutine init_itv_140

  real(kind(1.d0)) function itv_140()
    use tfcoil_variables, only: dr_tf_wp
    implicit none
    itv_140 = dr_tf_wp
  end function itv_140

  subroutine set_itv_140(ratio)
    use tfcoil_variables, only: dr_tf_wp
    implicit none
    real(kind(1.d0)) :: ratio
    dr_tf_wp = ratio
  end subroutine set_itv_140

  !---------------------------------

  subroutine init_itv_141
    !! <LI> (141) fcqt : TF coil quench temperature < tmax_croco (f-value for equation 74)
    use numerics, only: lablxc, boundl, boundu
    implicit none
    lablxc(141) = 'fcqt          '
    boundl(141) = 0.001D0
    boundu(141) = 1.000D0
  end subroutine init_itv_141

  real(kind(1.d0)) function itv_141()
    use constraint_variables, only: fcqt
    implicit none
    itv_141 = fcqt
  end function itv_141

  subroutine set_itv_141(ratio)
    use constraint_variables, only: fcqt
    implicit none
    real(kind(1.d0)) :: ratio
    fcqt = ratio
  end subroutine set_itv_141

  !---------------------------------

  subroutine init_itv_142
    !! <LI> (142) nesep : electron density at separatrix [m-3]
    use numerics, only: lablxc, boundl, boundu
    implicit none
    lablxc(142) = 'nesep         '
    boundl(142) = 1.00D17
    boundu(142) = 1.00D20
  end subroutine init_itv_142

  real(kind(1.d0)) function itv_142()
    use physics_variables, only: nesep
    implicit none
    itv_142 = nesep
  end function itv_142

  subroutine set_itv_142(ratio)
    use physics_variables, only: nesep
    implicit none
    real(kind(1.d0)) :: ratio
    nesep = ratio
  end subroutine set_itv_142

  !---------------------------------

  subroutine init_itv_143
    !! <LI> (143) f_copperA_m2 : TF coil current / copper area < Maximum value
    !!            (f-value for equation 75)
    use numerics, only: lablxc, boundl, boundu
    implicit none
    lablxc(143) = 'f_coppera_m2  '
    boundl(143) = 0.001D0
    boundu(143) = 1.000D0
  end subroutine init_itv_143

  real(kind(1.d0)) function itv_143()
    use rebco_variables, only: f_coppera_m2
    implicit none
    itv_143 = f_copperA_m2
  end function itv_143

  subroutine set_itv_143(ratio)
    use rebco_variables, only: f_coppera_m2
    implicit none
    real(kind(1.d0)) :: ratio
    f_coppera_m2 = ratio
  end subroutine set_itv_143

  !---------------------------------

  subroutine init_itv_144
    !! <LI> (144) fnesep : Eich critical electron density at separatrix
    use numerics, only: lablxc, boundl, boundu
    implicit none
    !!                     (f-value for constraint equation 76) [m-3]
    lablxc(144) = 'fnesep        '
    boundl(144) = 0.001D0
    boundu(144) = 1.000D0
  end subroutine init_itv_144

  real(kind(1.d0)) function itv_144()
    use constraint_variables, only: fnesep
    implicit none
    itv_144 = fnesep
  end function itv_144

  subroutine set_itv_144(ratio)
    use constraint_variables, only: fnesep
    implicit none
    real(kind(1.d0)) :: ratio
    fnesep = ratio
  end subroutine set_itv_144

  !---------------------------------

  subroutine init_itv_145
    !! <LI> (145) fgwped :  fraction of Greenwald density to set as pedestal-top density
    use numerics, only: lablxc, boundl, boundu
    implicit none
    lablxc(145) = 'fgwped        '
    boundl(145) = 0.500D0
    boundu(145) = 1.000D0
  end subroutine init_itv_145

  real(kind(1.d0)) function itv_145()
    use physics_variables, only: fgwped
    implicit none
    itv_145 = fgwped
  end function itv_145

  subroutine set_itv_145(ratio)
    use physics_variables, only: fgwped
    implicit none
    real(kind(1.d0)) :: ratio
    fgwped = ratio
  end subroutine set_itv_145

  !---------------------------------

  subroutine init_itv_146
    !! <LI> (146) fcpttf : F-value for TF coil current per turn limit (constraint equation 77)
    use numerics, only: lablxc, boundl, boundu
    implicit none
    lablxc(146) = 'fcpttf        '
    boundl(146) = 0.001D0
    boundu(146) = 1.000D0
  end subroutine init_itv_146

  real(kind(1.d0)) function itv_146()
    use constraint_variables, only: fcpttf
    implicit none
    itv_146 = fcpttf
  end function itv_146

  subroutine set_itv_146(ratio)
    use constraint_variables, only: fcpttf
    implicit none
    real(kind(1.d0)) :: ratio
    fcpttf = ratio
  end subroutine set_itv_146

  !---------------------------------

  subroutine init_itv_147
    !! <LI> (147) freinke : F-value for Reinke detachment criterion (constraint equation 78)
    use numerics, only: lablxc, boundl, boundu
    implicit none
    lablxc(147) = 'freinke       '
    boundl(147) = 0.001D0
    boundu(147) = 1.000D0
  end subroutine init_itv_147

  real(kind(1.d0)) function itv_147()
    use constraint_variables, only: freinke
    implicit none
    itv_147 = freinke
  end function itv_147

  subroutine set_itv_147(ratio)
    use constraint_variables, only: freinke
    implicit none
    real(kind(1.d0)) :: ratio
    freinke = ratio
  end subroutine set_itv_147

  !---------------------------------

  subroutine init_itv_148
    !! <LI> (148) fzactual : fraction of impurity at SOL with Reinke detachment criterion
    use numerics, only: lablxc, boundl, boundu
    implicit none
    lablxc(148) = 'fzactual      '
    boundl(148) = 1.00D-8
    boundu(148) = 1.000D0
  end subroutine init_itv_148

  real(kind(1.d0)) function itv_148()
    ! From Kallenbach model, should be reserved if the model is going to be added back
  end function itv_148

  subroutine set_itv_148(ratio)
    implicit none
    real(kind(1.d0)) :: ratio

    ! From Kallenbach model, should be reserved if the model is going to be added back
  end subroutine set_itv_148

  !---------------------------------

  subroutine init_itv_149
    !! <LI> (149) fbmaxcs : F-value for max peak CS field (con. 79, itvar 149)
    use numerics, only: lablxc, boundl, boundu
    implicit none
    lablxc(149) = 'fbmaxcs       '
    boundl(149) = 0.001D0
    boundu(149) = 1.000D0
  end subroutine init_itv_149

  real(kind(1.d0)) function itv_149()
    use pfcoil_variables, only: fbmaxcs
    implicit none
    itv_149 = fbmaxcs
  end function itv_149

  subroutine set_itv_149(ratio)
    use pfcoil_variables, only: fbmaxcs
    implicit none
    real(kind(1.d0)) :: ratio
    fbmaxcs = ratio
  end subroutine set_itv_149

  !---------------------------------

  real(kind(1.d0)) function itv_150()
    implicit none
    write(*,*) 'Iteration variable 150 no longer in use'
  end function itv_150

  !---------------------------------

  real(kind(1.d0)) function itv_151()
    implicit none
    write(*,*) 'Iteration variable 150 no longer in use'
  end function itv_151

  !---------------------------------

  subroutine init_itv_152
    !! <LI> (152) fbmaxcs : Ratio of separatrix density to Greenwald density
    use numerics, only: lablxc, boundl, boundu
    implicit none
    lablxc(152) = 'fgwsep        '
    boundl(152) = 0.001D0
    boundu(152) = 1.000D0
  end subroutine init_itv_152

  real(kind(1.d0)) function itv_152()
    use physics_variables, only: fgwsep
    implicit none
    itv_152 = fgwsep
  end function itv_152

  subroutine set_itv_152(ratio)
    use physics_variables, only: fgwsep
    implicit none
    real(kind(1.d0)) :: ratio
    fgwsep = ratio
  end subroutine set_itv_152

  !---------------------------------

  subroutine init_itv_153
    !! <LI> (153) fpdivlim : F-value for minimum pdivt (con. 80)
    use numerics, only: lablxc, boundl, boundu
    implicit none
    lablxc(153) = 'fpdivlim      '
    boundl(153) = 0.001D0
    boundu(153) = 1.000D0
  end subroutine init_itv_153

  real(kind(1.d0)) function itv_153()
    use physics_variables, only: fpdivlim
    implicit none
    itv_153 = fpdivlim
  end function itv_153

  subroutine set_itv_153(ratio)
    use physics_variables, only: fpdivlim
    implicit none
    real(kind(1.d0)) :: ratio
    fpdivlim = ratio
  end subroutine set_itv_153

  !---------------------------------

  subroutine init_itv_154
    !! <LI> (154) fne0 : F-value for ne(0) > ne(ped) (con. 81)
    use numerics, only: lablxc, boundl, boundu
    implicit none
    lablxc(154) = 'fne0          '
    boundl(154) = 0.001D0
    boundu(154) = 1.000D0
  end subroutine init_itv_154

  real(kind(1.d0)) function itv_154()
    use physics_variables, only: fne0
    implicit none
    itv_154 = fne0
  end function itv_154

  subroutine set_itv_154(ratio)
    use physics_variables, only: fne0
    implicit none
    real(kind(1.d0)) :: ratio
    fne0 = ratio
  end subroutine set_itv_154

  !---------------------------------

  subroutine init_itv_155
    !! <LI> (155) pfusife : IFE input fusion power (MW) (ifedrv=3 only)
    use numerics, only: lablxc, boundl, boundu
    implicit none
    lablxc(155) = 'pfusife       '
    boundl(155) = 5.000d2
    boundu(155) = 3.000d3
  end subroutine init_itv_155

  real(kind(1.d0)) function itv_155()
    use ife_variables, only: pfusife
    implicit none
    itv_155 = pfusife
  end function itv_155

  subroutine set_itv_155(ratio)
    use ife_variables, only: pfusife
    implicit none
    real(kind(1.d0)) :: ratio
    pfusife = ratio
  end subroutine set_itv_155

  !---------------------------------

  subroutine init_itv_156
    !! <LI> (156) rrin : Input IFE repetition rate (Hz) (ifedrv=3 only)
    use numerics, only: lablxc, boundl, boundu
    implicit none
    lablxc(156) = 'rrin          '
    boundl(156) = 1.000d0
    boundu(156) = 1.000d1
  end subroutine init_itv_156

  real(kind(1.d0)) function itv_156()
    use ife_variables, only: rrin
    implicit none
    itv_156 = rrin
  end function itv_156

  subroutine set_itv_156(ratio)
    use ife_variables, only: rrin
    implicit none
    real(kind(1.d0)) :: ratio
    rrin = ratio
  end subroutine set_itv_156

  !---------------------------------

  subroutine init_itv_157
    !! <LI> (157) fvssu : F-value for available to required start up flux (con. 51)
    use numerics, only: lablxc, boundl, boundu
    implicit none
    lablxc(157) = 'fvssu         '
    boundl(157) = 1.00d-3
    boundu(157) = 1.000d1
  end subroutine init_itv_157

  real(kind(1.d0)) function itv_157()
    use pfcoil_variables, only: fvssu
    implicit none
    itv_157 = fvssu
  end function itv_157

  subroutine set_itv_157(ratio)
    use pfcoil_variables, only: fvssu
    implicit none
    real(kind(1.d0)) :: ratio
    fvssu = ratio
  end subroutine set_itv_157

  !---------------------------------

  subroutine init_itv_158
    !! <LI> (158) croco_thick : Thickness of CroCo copper tube (m)
    use numerics, only: lablxc, boundl, boundu
    implicit none
    lablxc(158) = 'croco_thick   '
    boundl(158) = 1.0d-3
    boundu(158) = 1.0d-1
  end subroutine init_itv_158

  real(kind(1.d0)) function itv_158()
    use rebco_variables, only: croco_thick
    implicit none
    itv_158 = croco_thick
  end function itv_158

  subroutine set_itv_158(ratio)
    use rebco_variables, only: croco_thick
    implicit none
    real(kind(1.d0)) :: ratio
    croco_thick = ratio
  end subroutine set_itv_158

  !---------------------------------

  subroutine init_itv_159
    !! <LI> (159) ftoroidalgap : F-value for toroidalgap >  tftort constraint (con. 82)
    use numerics, only: lablxc, boundl, boundu
    lablxc(159) = 'ftoroidalgap  '
    boundl(159) = 1.0D-4
    boundu(159) = 1.0D0
  end subroutine init_itv_159

  real(kind(1.d0)) function itv_159()
    use tfcoil_variables, only: ftoroidalgap
    implicit none
    itv_159 = ftoroidalgap
  end function itv_159

  subroutine set_itv_159(ratio)
    use tfcoil_variables, only: ftoroidalgap
    implicit none
    real(kind(1.d0)) :: ratio
    ftoroidalgap = ratio
  end subroutine set_itv_159

  !---------------------------------

  subroutine init_itv_160
    !! <LI> (160) f_avspace (f-value for equation 83)
    use numerics, only: lablxc, boundu, boundl
    lablxc(160) = 'f_avspace     '
    boundl(160) = 0.010D0
    boundu(160) = 1.000D0
  end subroutine init_itv_160

  real(kind(1.d0)) function itv_160()
    use build_variables, only: f_avspace
    implicit none
    itv_160 = f_avspace
  end function itv_160

  subroutine set_itv_160(ratio)
    use build_variables, only: f_avspace
    implicit none
    real(kind(1.d0)) :: ratio
    f_avspace = ratio
  end subroutine set_itv_160

  !---------------------------------

  subroutine init_itv_161
    !! <LI> (161) fbetatry_lower (f-value for equation 84)
    use constraint_variables, only: fbetatry_lower
    use numerics, only: lablxc, boundl, boundu
    implicit none
    lablxc(161) = 'fbetatry_lower     '
    boundl(161) = 0.010D0
    boundu(161) = 1.000D0
  end subroutine init_itv_161

  real(kind(1.d0)) function itv_161()
    use constraint_variables, only: fbetatry_lower
    implicit none
    itv_161 = fbetatry_lower
  end function itv_161

  subroutine set_itv_161(ratio)
    use constraint_variables, only: fbetatry_lower
    real(kind(1.d0)) :: ratio
    fbetatry_lower = ratio
  end subroutine set_itv_161

  !---------------------------------

  subroutine init_itv_162
    !! <LI> (162) r_cp_top : Top outer radius of the centropost (ST only) (m)
    use numerics, only: lablxc, boundl, boundu
    implicit none
    lablxc(162) = 'r_cp_top           '
    boundl(162) = 0.0010D0
    boundu(162) = 10.000D0
  end subroutine init_itv_162

  real(kind(1.d0)) function itv_162()
  use build_variables, only: r_cp_top
    implicit none
    itv_162 = r_cp_top
  end function itv_162

  subroutine set_itv_162(ratio)
    use build_variables, only: r_cp_top
    real(kind(1.d0)) :: ratio
    r_cp_top = ratio
  end subroutine set_itv_162

  !---------------------------------

  subroutine init_itv_163
    !! <LI> (163)f_t_turn_tf : Top outer radius of the centropost (ST only) (m)
    use numerics, only: lablxc, boundl, boundu
    implicit none
    lablxc(163) = 'f_t_turn_tf        '
    boundl(163) = 0.0010D0
    boundu(163) = 1000.0D0
  end subroutine init_itv_163

  real(kind(1.d0)) function itv_163()
  use tfcoil_variables, only: f_t_turn_tf
    implicit none
    itv_163 = f_t_turn_tf
  end function itv_163

  subroutine set_itv_163(ratio)
    use tfcoil_variables, only: f_t_turn_tf
    real(kind(1.d0)) :: ratio
    f_t_turn_tf = ratio
  end subroutine set_itv_163

  !---------------------------------

  subroutine init_itv_164
    !! <LI> (164) f-value for maximum cryogenic plant power
    use numerics, only: lablxc, boundl, boundu
    implicit none
    lablxc(164) = 'f_crypmw         '
    boundl(164) = 0.001D0
    boundu(164) = 1.000D0
  end subroutine init_itv_164

  real(kind(1.d0)) function itv_164()
    use heat_transport_variables, only: f_crypmw
    implicit none
    itv_164 = f_crypmw
  end function itv_164

  subroutine set_itv_164(ratio)
    use heat_transport_variables, only: f_crypmw
    implicit none
    real(kind(1.d0)) :: ratio
    f_crypmw = ratio
  end subroutine set_itv_164

  !---------------------------------

  subroutine init_itv_165
    !! <LI> (165) f-value for maximum TF coil strain
    use numerics, only: lablxc, boundl, boundu
    implicit none
    lablxc(165) = 'fstr_wp       '
    boundl(165) = 1.0d-9
    boundu(165) = 1.0d0
  end subroutine init_itv_165

  real(kind(1.d0)) function itv_165()
    use constraint_variables, only: fstr_wp
    implicit none
    itv_165 = fstr_wp
  end function itv_165

  subroutine set_itv_165(ratio)
    use constraint_variables, only: fstr_wp
    implicit none
    real(kind(1.d0)) :: ratio
    fstr_wp = ratio
  end subroutine  set_itv_165

  !---------------------------------

  subroutine init_itv_166
    !! <LI> (166) f_copperaoh_m2 : CS coil current / copper area < Maximum value
    !! (f-value for equation 89)
    use numerics, only: lablxc, boundl, boundu
    implicit none
    lablxc(166) = 'f_copperaoh_m2'
    boundl(166) = 0.001D0
    boundu(166) = 1.000D0
  end subroutine init_itv_166

  real(kind(1.d0)) function itv_166()
  use rebco_variables, only: f_copperaoh_m2
    implicit none
    itv_166 = f_copperaoh_m2
  end function itv_166

  subroutine set_itv_166(ratio)
    use rebco_variables, only: f_copperaoh_m2
    implicit none
    real(kind(1.d0)) :: ratio
    f_copperaoh_m2 = ratio
  end subroutine  set_itv_166

  !---------------------------------

  subroutine init_itv_167
    !! <LI> (167) fncycle : f-value for minimum CS coil stress load cycles
    use numerics, only: lablxc, boundl, boundu
    implicit none
    lablxc(167) = 'fncycle       '
    boundl(167) = 1.0d-8
    boundu(167) = 1.0d0
  end subroutine init_itv_167

  real(kind(1.d0)) function itv_167()
    use constraint_variables, only: fncycle
    implicit none
    itv_167 = fncycle
  end function itv_167

  subroutine set_itv_167(ratio)
    use constraint_variables, only: fncycle
    implicit none
    real(kind(1.d0)) :: ratio
    fncycle = ratio
  end subroutine  set_itv_167

  !---------------------------------

  subroutine init_itv_168
    !! <LI> (168) fecrh_ignition (f-value for equation 91)
    use numerics, only: lablxc, boundl, boundu
    implicit none
    lablxc(168) = 'fecrh_ignition     '
    boundl(168) = 0.010D0
    boundu(168) = 2.000D0
  end subroutine init_itv_168

  real(kind(1.d0)) function itv_168()
    use constraint_variables, only: fecrh_ignition
    implicit none
    itv_168 = fecrh_ignition
  end function itv_168

  subroutine set_itv_168(ratio)
    use constraint_variables, only: fecrh_ignition
    implicit none
    real(kind(1.d0)) :: ratio
    fecrh_ignition = ratio
  end subroutine set_itv_168

  !---------------------------------

  subroutine init_itv_169
    !! <LI> (169) te0_ecrh_achievable (iteration parameter for equation 91)
    use numerics, only: lablxc, boundl, boundu
    implicit none
    lablxc(169) = 'te0_ecrh_achievable'
    boundl(169) = 1.0d0
    boundu(169) = 1.0d3
  end subroutine init_itv_169

  real(kind(1.d0)) function itv_169()
    use stellarator_variables, only: te0_ecrh_achievable
    implicit none
    itv_169 = te0_ecrh_achievable
  end function itv_169

  subroutine set_itv_169(ratio)
    use stellarator_variables, only: te0_ecrh_achievable
    implicit none
    real(kind(1.d0)) :: ratio
    te0_ecrh_achievable = ratio
  end subroutine set_itv_169

  subroutine init_itv_170
    !! <LI> (170) beta_div : field line angle wrt divertor target plate (degrees)
    use numerics, only: lablxc, boundl, boundu
    implicit none
    lablxc(170) = 'beta_div      '
    boundl(170) = 0.49
    boundu(170) = 5.01
  end subroutine init_itv_170

  real(kind(1.d0)) function itv_170()
    use divertor_variables, only: beta_div
    implicit none
    itv_170 = beta_div
  end function itv_170

  subroutine set_itv_170(ratio)
    use divertor_variables, only: beta_div
    implicit none
    real(kind(1.d0)) :: ratio
    beta_div = ratio
  end subroutine set_itv_170

  !---------------------------------
  ! DUMMY variables below here
  !---------------------------------

  subroutine init_itv_171
    !! <LI> (171) DUMMY : Description
    use numerics, only: lablxc, boundl, boundu
    implicit none
    lablxc(171) = 'DUMMY         '
    boundl(171) = 1.0d-99
    boundu(171) = 1.0d99
  end subroutine init_itv_171

  real(kind(1.d0)) function itv_171()
    implicit none
    itv_171 = DUMMY
  end function itv_171

  subroutine set_itv_171(ratio)
    implicit none
    real(kind(1.d0)) :: ratio
    DUMMY = ratio
  end subroutine set_itv_171

  !---------------------------------

  subroutine init_itv_172
    !! <LI> (172) DUMMY : Description
    use numerics, only: lablxc, boundl, boundu
    implicit none
    lablxc(172) = 'DUMMY         '
    boundl(172) = 1.0d-99
    boundu(172) = 1.0d99
  end subroutine init_itv_172

  real(kind(1.d0)) function itv_172()
    implicit none
    itv_172 = DUMMY
  end function itv_172

  subroutine set_itv_172(ratio)
    implicit none
    real(kind(1.d0)) :: ratio
    DUMMY = ratio
  end subroutine set_itv_172

  !---------------------------------

  subroutine init_itv_173
    !! <LI> (173) DUMMY : Description
    use numerics, only: lablxc, boundl, boundu
    implicit none
    lablxc(173) = 'DUMMY         '
    boundl(173) = 1.0d-99
    boundu(173) = 1.0d99
  end subroutine init_itv_173

  real(kind(1.d0)) function itv_173()
    implicit none
    itv_173 = DUMMY
  end function itv_173

  subroutine set_itv_173(ratio)
    implicit none
    real(kind(1.d0)) :: ratio
    DUMMY = ratio
  end subroutine set_itv_173

  !---------------------------------

  subroutine init_itv_174
    !! <LI> (174) DUMMY : Description
    use numerics, only: lablxc, boundl, boundu
    implicit none
    lablxc(174) = 'DUMMY         '
    boundl(174) = 1.0d-99
    boundu(174) = 1.0d99
  end subroutine init_itv_174

  real(kind(1.d0)) function itv_174()
    implicit none
    itv_174 = DUMMY
  end function itv_174

  subroutine set_itv_174(ratio)
    implicit none
    real(kind(1.d0)) :: ratio
    DUMMY = ratio
  end subroutine set_itv_174

  !---------------------------------

  subroutine init_itv_175
    !! <LI> (175) DUMMY : Description
    use numerics, only: lablxc, boundl, boundu
    implicit none
    lablxc(175) = 'DUMMY         '
    boundl(175) = 1.0d-99
    boundu(175) = 1.0d99
  end subroutine init_itv_175

  real(kind(1.d0)) function itv_175()
    implicit none
    itv_175 = DUMMY
  end function itv_175

  subroutine set_itv_175(ratio)
    implicit none
    real(kind(1.d0)) :: ratio
    DUMMY = ratio
  end subroutine set_itv_175

  !------------------------------------------------------

  subroutine loadxc
    !! Routine to load the physics and engineering variables into the
    !! optimisation variables array
    !! author: P J Knight, CCFE, Culham Science Centre
    !! author: J Morris, CCFE, Culham Science Centre
    !! None
    !! This subroutine loads the physics and engineering variables
    !! into the optimisation variables array <CODE>XCM</CODE>.
    !
    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    use constants, only: nout
    use maths_library, only: variable_error
    use error_handling, only: idiags, fdiags, report_error
    use numerics, only: nvar, xcm, ixc, name_xc, lablxc, scafc, scale
    use physics_variables, only: icurr
    use global_variables, only: vlabel
    implicit none

    !  Local variables
    integer :: i
    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    do i = 1,nvar

       select case (ixc(i))
           case (1);  xcm(i) = itv_1()
           case (2);  xcm(i) = itv_2()
           case (3);  xcm(i) = itv_3()
           case (4);  xcm(i) = itv_4()
           case (5);  xcm(i) = itv_5()
           case (6);  xcm(i) = itv_6()
           case (7);  xcm(i) = itv_7()
           case (8);  xcm(i) = itv_8()
           case (9);  xcm(i) = itv_9()
           case (10);  xcm(i) = itv_10()
           case (11);  xcm(i) = itv_11()
           case (12);  xcm(i) = itv_12()
           case (13);  xcm(i) = itv_13()
           case (14);  xcm(i) = itv_14()
           case (15);  xcm(i) = itv_15()
           case (16);  xcm(i) = itv_16()
           case (17);  xcm(i) = itv_17()
           case (18);  xcm(i) = itv_18()
           case (19);  xcm(i) = itv_19()
           case (20);  xcm(i) = itv_20()
           case (21);  xcm(i) = itv_21()
           case (22);  xcm(i) = itv_22()
           case (23);  xcm(i) = itv_23()
           case (24);  xcm(i) = itv_24()
           case (25);  xcm(i) = itv_25()
           case (26);  xcm(i) = itv_26()
           case (27);  xcm(i) = itv_27()
           case (28);  xcm(i) = itv_28()
           case (29);  xcm(i) = itv_29()
           case (30);  xcm(i) = itv_30()
           case (31);  xcm(i) = itv_31()
           case (32);  xcm(i) = itv_32()
           case (33);  xcm(i) = itv_33()
           case (34);  xcm(i) = itv_34()
           case (35);  xcm(i) = itv_35()
           case (36);  xcm(i) = itv_36()
           case (37);  xcm(i) = itv_37()
           case (38);  xcm(i) = itv_38()
           case (39);  xcm(i) = itv_39()
           case (40);  xcm(i) = itv_40()
           case (41);  xcm(i) = itv_41()
           case (42);  xcm(i) = itv_42()
           case (43);
           case (44);  xcm(i) = itv_44()
           case (45);  xcm(i) = itv_45()
           case (46);  xcm(i) = itv_46()
           case (47);  xcm(i) = itv_47()
           case (48);  xcm(i) = itv_48()
           case (49);  xcm(i) = itv_49()
           case (50);  xcm(i) = itv_50()
           case (51);  xcm(i) = itv_51()
           case (52);  xcm(i) = itv_52()
           case (53);  xcm(i) = itv_53()
           case (54);  xcm(i) = itv_54()
           case (55);  xcm(i) = itv_55()
           case (56);  xcm(i) = itv_56()
           case (57);  xcm(i) = itv_57()
           case (58);  xcm(i) = itv_58()
           case (59);  xcm(i) = itv_59()
           case (60);  xcm(i) = itv_60()
           case (61);  xcm(i) = itv_61()
           case (62);  xcm(i) = itv_62()
           case (63);  xcm(i) = itv_63()
           case (64);  xcm(i) = itv_64()
           case (65);  xcm(i) = itv_65()
           case (66);  xcm(i) = itv_66()
           case (67);  xcm(i) = itv_67()
           case (68);  xcm(i) = itv_68()
           case (69);  xcm(i) = itv_69()
           case (70);  xcm(i) = itv_70()
           case (71);  xcm(i) = itv_71()
           case (72);  xcm(i) = itv_72()
           case (73);  xcm(i) = itv_73()
           case (74);  xcm(i) = itv_74()
           case (75);  xcm(i) = itv_75()
           case (76);  xcm(i) = itv_76()
           case (77);  xcm(i) = itv_77()
           case (78);  xcm(i) = itv_78()
           case (79);  xcm(i) = itv_79()
           case (80);  xcm(i) = itv_80()
           case (81);  xcm(i) = itv_81()
           case (82);  xcm(i) = itv_82()
           case (83);  xcm(i) = itv_83()
           case (84);  xcm(i) = itv_84()
           case (85);  xcm(i) = itv_85()
           case (86);  xcm(i) = itv_86()
           case (87);
           case (88);
           case (89);  xcm(i) = itv_89()
           case (90);  xcm(i) = itv_90()
           case (91);  xcm(i) = itv_91()
           case (92);  xcm(i) = itv_92()
           case (93);  xcm(i) = itv_93()
           case (94);  xcm(i) = itv_94()
           case (95);  xcm(i) = itv_95()
           case (96);  xcm(i) = itv_96()
           case (97);  xcm(i) = itv_97()
           case (98);  xcm(i) = itv_98()
           case (99);
           case (100);
           case (101);
           case (102);
           case (103);  xcm(i) = itv_103()
           case (104);  xcm(i) = itv_104()
           case (105);  xcm(i) = itv_105()
           case (106);  xcm(i) = itv_106()
           case (107);  xcm(i) = itv_107()
           case (108);  xcm(i) = itv_108()
           case (109);  xcm(i) = itv_109()
           case (110);  xcm(i) = itv_110()
           case (111);  xcm(i) = itv_111()
           case (112);  xcm(i) = itv_112()
           case (113);  xcm(i) = itv_113()
           case (114);  xcm(i) = itv_114()
           case (115);  xcm(i) = itv_115()
           case (116);  xcm(i) = itv_116()
           case (117);  xcm(i) = itv_117()
           case (118);  xcm(i) = itv_118()
           case (119);  xcm(i) = itv_119()
           case (120);  xcm(i) = itv_120()
           case (121);  xcm(i) = itv_121()
           case (122);  xcm(i) = itv_122()
           case (123);  xcm(i) = itv_123()
           case (124);  xcm(i) = itv_124()
           case (125);  xcm(i) = itv_125()
           case (126);  xcm(i) = itv_126()
           case (127);  xcm(i) = itv_127()
           case (128);  xcm(i) = itv_128()
           case (129);  xcm(i) = itv_129()
           case (130);  xcm(i) = itv_130()
           case (131);  xcm(i) = itv_131()
           case (132);  xcm(i) = itv_132()
           case (133);  xcm(i) = itv_133()
           case (134);  xcm(i) = itv_134()
           case (135);  xcm(i) = itv_135()
           case (136);  xcm(i) = itv_136()
           case (137);  xcm(i) = itv_137()
           case (138);  xcm(i) = itv_138()
           case (139);  xcm(i) = itv_139()
           case (140);  xcm(i) = itv_140()
           case (141);  xcm(i) = itv_141()
           case (142);  xcm(i) = itv_142()
           case (143);  xcm(i) = itv_143()
           case (144);  xcm(i) = itv_144()
           case (145);  xcm(i) = itv_145()
           case (146);  xcm(i) = itv_146()
           case (147);  xcm(i) = itv_147()
           case (148);  xcm(i) = itv_148()
           case (149);  xcm(i) = itv_149()
           case (150);
           case (151);
           case (152);  xcm(i) = itv_152()
           case (153);  xcm(i) = itv_153()
           case (154);  xcm(i) = itv_154()
           case (155);  xcm(i) = itv_155()
           case (156);  xcm(i) = itv_156()
           case (157);  xcm(i) = itv_157()
           case (158);  xcm(i) = itv_158()
           case (159);  xcm(i) = itv_159()
           case (160);  xcm(i) = itv_160()
           case (161);  xcm(i) = itv_161()
           case (162);  xcm(i) = itv_162()
           case (163);  xcm(i) = itv_163()
           case (164);  xcm(i) = itv_164()
           case (165);  xcm(i) = itv_165()
           case (166);  xcm(i) = itv_166()
           case (167);  xcm(i) = itv_167()
           case (168);  xcm(i) = itv_168()
           case (169);  xcm(i) = itv_169()
            ! DUMMY Cases
           case (170);  xcm(i) = itv_170()
           case (171);  xcm(i) = itv_171()
           case (172);  xcm(i) = itv_172()
           case (173);  xcm(i) = itv_173()
           case (174);  xcm(i) = itv_174()
           case (175);  xcm(i) = itv_175()

       case default
          idiags(1) = i ; idiags(2) = ixc(i)
          call report_error(54)

       end select


        ! Simple list of iteration variable names
        name_xc(i) = lablxc(ixc(i))
        ! Note that iteration variable 18 has more than one name:
        if ((ixc(i) == 18).and.(icurr /= 2)) name_xc(i) = 'q95'
        if ((ixc(i) == 18).and.(icurr == 2)) name_xc(i) = 'qbar'


         ! MDK Check if sweep variable is also an iteration variable
         if (name_xc(i) == vlabel) then
              write(nout,*) 'WARNING: The sweep variable is also an iteration variable.'
              write(nout,*) 'The values of the sweep variable will be overwritten by the optimiser.'
              write(*,*) 'WARNING: The sweep variable is also an iteration variable.'
         end if

       !  Check that no iteration variable is zero

       if (abs(xcm(i)) <= 1.0D-12) then
          idiags(1) = i ; idiags(2) = ixc(i)
          write(*,*) 'Iteration variable ',ixc(i),'(',trim(lablxc(ixc(i))),') is zero.'
          call report_error(55)
       end if

       !  Crude method of catching NaN errors

       !if ( (abs(xcm(i)) > 9.99D99).or.(xcm(i) /= xcm(i)) ) then
       if ( variable_error(xcm(i)) ) then
          idiags(1) = i ; idiags(2) = ixc(i) ; fdiags(1) = xcm(i)
          call report_error(56)
       end if

    end do

    do i = 1,nvar
       if (abs(xcm(i)) > 1.0D-99) then
          scale(i) = 1.0D0/xcm(i)
       else
          scale(i) = 1.0D0
       end if
       scafc(i) = 1.0D0/scale(i)
       xcm(i) = xcm(i)*scale(i)
    end do

  end subroutine loadxc
  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine convxc(xc,nn)

    !! Routine to convert scaled iteration variables back to
    !! their real values
    !! author: P J Knight, CCFE, Culham Science Centre
    !! author: J Morris, CCFE, Culham Science Centre
    !! xc(ipnvars) : input/output real array : scaled iteration variable values
    !! nn : input integer : number of iteration variables
    !! This subroutine converts the scaled iteration variables back to
    !! their real values.
    !
    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    use error_handling, only: idiags, fdiags, report_error
    use numerics, only: ipnvars, scale, ixc, lablxc
    use maths_library, only: variable_error
#ifndef dp
    use, intrinsic :: iso_fortran_env, only: dp=>real64
#endif

    implicit none

    !  Arguments

    integer, intent(in) :: nn
    real(dp), dimension(:), intent(in) :: xc

    !  Local variables

    integer :: i
    real(dp)::ratio

    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    do i = 1,nn
       ratio = xc(i)/scale(i)

       select case (ixc(i))
           case (1);  call set_itv_1(ratio)
           case (2);  call set_itv_2(ratio)
           case (3);  call set_itv_3(ratio)
           case (4);  call set_itv_4(ratio)
           case (5);  call set_itv_5(ratio)
           case (6);  call set_itv_6(ratio)
           case (7);  call set_itv_7(ratio)
           case (8);  call set_itv_8(ratio)
           case (9);  call set_itv_9(ratio)
           case (10);  call set_itv_10(ratio)
           case (11);  call set_itv_11(ratio)
           case (12);  call set_itv_12(ratio)
           case (13);  call set_itv_13(ratio)
           case (14);  call set_itv_14(ratio)
           case (15);  call set_itv_15(ratio)
           case (16);  call set_itv_16(ratio)
           case (17);  call set_itv_17(ratio)
           case (18);  call set_itv_18(ratio)
           case (19);  call set_itv_19(ratio)
           case (20);  call set_itv_20(ratio)
           case (21);  call set_itv_21(ratio)
           case (22);
           case (23);  call set_itv_23(ratio)
           case (24);
           case (25);  call set_itv_25(ratio)
           case (26);  call set_itv_26(ratio)
           case (27);  call set_itv_27(ratio)
           case (28);  call set_itv_28(ratio)
           case (29);  call set_itv_29(ratio)
           case (30);  call set_itv_30(ratio)
           case (31);  call set_itv_31(ratio)
           case (32);  call set_itv_32(ratio)
           case (33);  call set_itv_33(ratio)
           case (34);  call set_itv_34(ratio)
           case (35);  call set_itv_35(ratio)
           case (36);  call set_itv_36(ratio)
           case (37);  call set_itv_37(ratio)
           case (38);  call set_itv_38(ratio)
           case (39);  call set_itv_39(ratio)
           case (40);  call set_itv_40(ratio)
           case (41);  call set_itv_41(ratio)
           case (42);  call set_itv_42(ratio)
           case (43);
           case (44);  call set_itv_44(ratio)
           case (45);  call set_itv_45(ratio)
           case (46);  call set_itv_46(ratio)
           case (47);  call set_itv_47(ratio)
           case (48);  call set_itv_48(ratio)
           case (49);  call set_itv_49(ratio)
           case (50);  call set_itv_50(ratio)
           case (51);  call set_itv_51(ratio)
           case (52);  call set_itv_52(ratio)
           case (53);  call set_itv_53(ratio)
           case (54);  call set_itv_54(ratio)
           case (55);
           case (56);  call set_itv_56(ratio)
           case (57);  call set_itv_57(ratio)
           case (58);  call set_itv_58(ratio)
           case (59);  call set_itv_59(ratio)
           case (60);  call set_itv_60(ratio)
           case (61);  call set_itv_61(ratio)
           case (62);  call set_itv_62(ratio)
           case (63);  call set_itv_63(ratio)
           case (64);  call set_itv_64(ratio)
           case (65);  call set_itv_65(ratio)
           case (66);  call set_itv_66(ratio)
           case (67);  call set_itv_67(ratio)
           case (68);  call set_itv_68(ratio)
           case (69);  call set_itv_69(ratio)
           case (70);  call set_itv_70(ratio)
           case (71);  call set_itv_71(ratio)
           case (72);  call set_itv_72(ratio)
           case (73);  call set_itv_73(ratio)
           case (74);  call set_itv_74(ratio)
           case (75);  call set_itv_75(ratio)
           case (76);
           case (77);
           case (78);
           case (79);  call set_itv_79(ratio)
           case (80);
           case (81);  call set_itv_81(ratio)
           case (82);  call set_itv_82(ratio)
           case (83);  call set_itv_83(ratio)
           case (84);  call set_itv_84(ratio)
           case (85);  call set_itv_85(ratio)
           case (86);  call set_itv_86(ratio)
           case (87);
           case (88);
           case (89);  call set_itv_89(ratio)
           case (90);  call set_itv_90(ratio)
           case (91);  call set_itv_91(ratio)
           case (92);  call set_itv_92(ratio)
           case (93);  call set_itv_93(ratio)
           case (94);  call set_itv_94(ratio)
           case (95);  call set_itv_95(ratio)
           case (96);  call set_itv_96(ratio)
           case (97);  call set_itv_97(ratio)
           case (98);  call set_itv_98(ratio)
           case (99);
           case (100);
           case (101);
           case (102);
           case (103);  call set_itv_103(ratio)
           case (104);  call set_itv_104(ratio)
           case (105);  call set_itv_105(ratio)
           case (106);  call set_itv_106(ratio)
           case (107);  call set_itv_107(ratio)
           case (108);  call set_itv_108(ratio)
           case (109);  call set_itv_109(ratio)
           case (110);  call set_itv_110(ratio)
           case (111);  call set_itv_111(ratio)
           case (112);  call set_itv_112(ratio)
           case (113);  call set_itv_113(ratio)
           case (114);  call set_itv_114(ratio)
           case (115);  call set_itv_115(ratio)
           case (116);  call set_itv_116(ratio)
           case (117);  call set_itv_117(ratio)
           case (118);  call set_itv_118(ratio)
           case (119);  call set_itv_119(ratio)
           case (120);  call set_itv_120(ratio)
           case (121);  call set_itv_121(ratio)
           case (122);  call set_itv_122(ratio)
           case (123);  call set_itv_123(ratio)
           case (124);  call set_itv_124(ratio)
           case (125);  call set_itv_125(ratio)
           case (126);  call set_itv_126(ratio)
           case (127);  call set_itv_127(ratio)
           case (128);  call set_itv_128(ratio)
           case (129);  call set_itv_129(ratio)
           case (130);  call set_itv_130(ratio)
           case (131);  call set_itv_131(ratio)
           case (132);  call set_itv_132(ratio)
           case (133);  call set_itv_133(ratio)
           case (134);  call set_itv_134(ratio)
           case (135);  call set_itv_135(ratio)
           case (136);  call set_itv_136(ratio)
           case (137);  call set_itv_137(ratio)
           case (138);  call set_itv_138(ratio)
           case (139);  call set_itv_139(ratio)
           case (140);  call set_itv_140(ratio)
           case (141);  call set_itv_141(ratio)
           case (142);  call set_itv_142(ratio)
           case (143);  call set_itv_143(ratio)
           case (144);  call set_itv_144(ratio)
           case (145);  call set_itv_145(ratio)
           case (146);  call set_itv_146(ratio)
           case (147);  call set_itv_147(ratio)
           case (148);  call set_itv_148(ratio)
           case (149);  call set_itv_149(ratio)
           case (150);
           case (151);
           case (152);  call set_itv_152(ratio)
           case (153);  call set_itv_153(ratio)
           case (154);  call set_itv_154(ratio)
           case (155);  call set_itv_155(ratio)
           case (156);  call set_itv_156(ratio)
           case (157);  call set_itv_157(ratio)
           case (158);  call set_itv_158(ratio)
           case (159);  call set_itv_159(ratio)
           case (160);  call set_itv_160(ratio)
           case (161);  call set_itv_161(ratio)
           case (162);  call set_itv_162(ratio)
           case (163);  call set_itv_163(ratio)
           case (164);  call set_itv_164(ratio)
           case (165);  call set_itv_165(ratio)
           case (166);  call set_itv_166(ratio)
           case (167);  call set_itv_167(ratio)
           case (168);  call set_itv_168(ratio)
           case (169);  call set_itv_169(ratio)
            ! DUMMY Cases
           case (170);  call set_itv_170(ratio)
           case (171);  call set_itv_171(ratio)
           case (172);  call set_itv_172(ratio)
           case (173);  call set_itv_173(ratio)
           case (174);  call set_itv_174(ratio)
           case (175);  call set_itv_175(ratio)

           case default
              call report_error(57)

       end select

       !  Check that no iteration variable is zero

       if (abs(xc(i)) <= 1.0D-12) then
          idiags(1) = i ; idiags(2) = ixc(i)
          write(*,*) 'Iteration variable ',ixc(i),'(',trim(lablxc(ixc(i))),') is zero.'
          call report_error(58)
       end if

       !  Crude method of catching NaN errors

       !if ((abs(xc(i)) > 9.99D99).or.(xc(i) /= xc(i))) then
       if(variable_error(xc(i)))then
          idiags(1) = i ; idiags(2) = ixc(i) ; fdiags(1) = xc(i)
          call report_error(59)
       end if

       if (abs(scale(i)) < 1.0D-99) then
          idiags(1) = i ; idiags(2) = ixc(i)
          call report_error(60)
       end if

    end do

  end subroutine convxc

  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine boundxc

    !! Routine to convert variable bounds to their real values
    !! author: P J Knight, CCFE, Culham Science Centre
    !! None
    !! This subroutine converts the scaled iteration variable bounds
    !! back to their real values.
    !! AEA FUS 251: A User's Guide to the PROCESS Systems Code
    !
    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    use numerics, only: nvar, bondl, bondu, scale, ixc, boundl, boundu

    implicit none

    !  Local variables

    integer :: i

    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    do i = 1,nvar
       bondl(i) = boundl(ixc(i))*scale(i)
       bondu(i) = boundu(ixc(i))*scale(i)
    end do

  end subroutine boundxc

  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
end module define_iteration_variables
! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

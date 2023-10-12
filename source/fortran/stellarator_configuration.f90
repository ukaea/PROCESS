module stellarator_configuration
    !! author: J Lion, IPP Greifswald
    !! Module containing defining parameters for a stellarator
    !!
    !! This module contains a set of constants that defines a
    !! stellarator configuration. These parameters are based on external
    !! calculations and are hardcoded right now into this module. There will be
    !! the possibiltiy to set them via an input file in the future.
    !! The list below will be modified in further commits.
    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#ifndef dp
   use, intrinsic :: iso_fortran_env, only: dp=>real64
#endif

    implicit none

    character (len = 20) :: stella_config_name
   ! Name of the configuration

   integer stella_config_symmetry
   !  Number of coils [1]

   integer stella_config_coilspermodule
   !  Coils per module [1]

   real(dp)  stella_config_rmajor_ref
   !  Reference Point for major radius where all the other variables are determined [m]

   real(dp)  stella_config_rminor_ref
   !  Reference Point for minor radius where all the other variables are determined [m]

   real(dp)  stella_config_coil_rmajor
   !  Reference Point for coil major radius [m]

   real(dp)  stella_config_coil_rminor
   !  Reference Point for coil minor radius [m]

   real(dp)  stella_config_aspect_ref
   !  Reference Point for aspect ratio where all the other variables are determined [1]

   real(dp)  stella_config_bt_ref
   !  Reference Point for toroidal b where all the other variables are determined [T]

   real(dp) stella_config_WP_area
   !  Winding pack area at the reference point [m^2]

   real(dp)  stella_config_WP_bmax
   !  The maximal magnetic field in the winding pack at the reference size of the winding pack [T]

   real(dp) stella_config_i0
   !  Coil current needed for b0 at the reference point [MA]

   real(dp) stella_config_a1
   !  Magnetic field fit parameter a1 (for the maximal field on the coils) [1]

   real(dp) stella_config_a2
   !  Magnetic field fit parameter a2 [1]

   real(dp) stella_config_dmin
   !  Minimal intercoil distance at the reference point [m]

   real(dp) stella_config_inductance
   !  inductance at the reference point [H]

   real(dp) stella_config_coilsurface
   !  Coil surface at the reference point [m2]

   real(dp) stella_config_coillength
   !  Total coil length at the reference point [m]

   real(dp) stella_config_max_portsize_width
   !  Port size in toroidal direction at the reference point [m]

   real(dp) stella_config_maximal_coil_height
   !  The maximal coil height at reference point. [m]

   real(dp) stella_config_min_plasma_coil_distance
   !  The minimal distance between coil and plasma at the reference point [m]

   real(dp) stella_config_derivative_min_LCFS_coils_dist
   !  The derivative of min_plasma_coil_distance wrt to the minor plasma radius at the reference point [1]

   real(dp) stella_config_plasma_volume
   !  The plasma volume at the reference point. Scales as a*R^2. [m^3]

   real(dp) stella_config_plasma_surface
   !  The plasma surface a the reference point. [m^2]

   real(dp) stella_config_WP_ratio
   !  Ratio radial to toroidal length of the winding pack. (a1 and a2 should be calculated using this value) [1]

   real(dp) stella_config_max_force_density
   !  Maximal toroidal and radially averaged force density at reference point in a WP cross section [MN/m^3]

   real(dp) stella_config_max_force_density_MNm
   !  Maximal integrated force density at reference point in a WP cross section [MN/m]

   real(dp) stella_config_min_bend_radius
   !  Minimal bending radius at reference point [m]

   real(dp) stella_config_epseff
   !  Maximal epsilon effective in the core region [1]

   real(dp) stella_config_max_lateral_force_density
   !  Maximal lateral force density of the coil set [MN/m]

   real(dp) stella_config_max_radial_force_density
   !  Maximal radial force density of the coil set [MN/m]

   real(dp) stella_config_centering_force_max_MN
   !  Maximal centering force of a coil in the coil set [MN]

   real(dp) stella_config_centering_force_min_MN
   !  Minimal centering force of a coil in the coil set (negative means pointing outwards) [MN]

   real(dp) stella_config_centering_force_avg_MN
   !  Average centering force the coils in the coil set [MN/coil]

   real(dp) :: stella_config_neutron_peakfactor
   !  The neutron peaking factor determined through inhomogeneities on the stellarator wall (qmax/qavg) [1]

   real(dp), dimension(:), allocatable :: sc_D11_star_mono_input
   !  The monoenergetic radial transport coefficients normalized by the plateau value.

   real(dp), dimension(:), allocatable :: sc_nu_star_mono_input
   !  The monoenergetic radial transport coefficients normalized by the plateau value.

 contains

    subroutine new_stella_config(index)
       integer, intent(in) :: index

       select case (index)

          ! This is the istell case switch:
          ! istell = 1: Helias5 machine
          ! istell = 2: Helias4 machine
          ! istell = 3: Helias3 machine
          ! istell = 4: w7x30 machine
          ! istell = 5: w7x50 machine
          ! istell = 6: Init from json

          ! All parameters set here are prelimnirary versions and might be changed in further commits

          case(1)
             ! Helias5 Machine
             ! The values are given at the reference point

             stella_config_name = "Helias 5b"

             stella_config_rmajor_ref = 22.2D0
             stella_config_rminor_ref = 1.80D0
             stella_config_aspect_ref = 12.33D0

             ! Coil radii
             stella_config_coil_rmajor = 22.44D0
             stella_config_coil_rminor = 4.76D0

             stella_config_bt_ref = 5.6D0
             stella_config_WP_area = 0.8d0*0.6d0
             stella_config_WP_bmax = 11.44d0

             stella_config_symmetry = 5
             stella_config_coilspermodule = 10

             stella_config_a1 = 0.688D0
             stella_config_a2 = 0.025D0

             stella_config_plasma_volume = 1422.63D0  ! This value is for Helias 5
             stella_config_dmin = 0.84D0
             stella_config_max_portsize_width = 2.12D0

             stella_config_plasma_surface = 1960.0D0 ! Plasma Surface

             stella_config_maximal_coil_height = 12.7D0 ! [m] Full height max point to min point

             stella_config_coilsurface = 4817.7D0 ! Coil surface, dimensionfull. At reference point

             stella_config_coillength = 1680.0D0 ! Central filament length of machine with outer radius 1m.

             stella_config_I0 = 13.06D0 ! Coil Current needed to produce 1T on axis in [MA] at outer radius 1m
             stella_config_inductance = 1655.76D-6 ! inductance in muH

             stella_config_WP_ratio = 1.2D0 ! The fit values in stellarator config class should be calculated using this value.

             stella_config_max_force_density = 120.0d0 ! [MN/m^3]
             stella_config_max_force_density_MNm = 98.0d0 ! [MN/m]

             stella_config_max_lateral_force_density = 92.4d0 ! [MN/m^3]
             stella_config_max_radial_force_density = 113.5d0   ! [MN/m^3]

             stella_config_centering_force_max_MN = 189.5d0
             stella_config_centering_force_min_MN = -55.7d0
             stella_config_centering_force_avg_MN = 93.0d0

             stella_config_min_plasma_coil_distance = 1.9d0
             stella_config_derivative_min_LCFS_coils_dist = -1.0d0 ! this is approximated for now

             stella_config_min_bend_radius = 1.0d0 ! [m]

             stella_config_neutron_peakfactor = 1.6d0

             stella_config_epseff = 0.015d0


             if (allocated(sc_D11_star_mono_input)) deallocate(sc_D11_star_mono_input)
             if (allocated(sc_nu_star_mono_input)) deallocate(sc_nu_star_mono_input)
             allocate(sc_D11_star_mono_input(10))
             allocate(sc_nu_star_mono_input(10))

             sc_D11_star_mono_input = (/1,1,1,1,1,1,1,1,1,1/)
             sc_nu_star_mono_input = (/1d-8,1d-7,1d-6,1d-5,1d-4,1d-3,1d-2,1d-1,1d0,1d1/)




          case(2)
             ! Helias4 Machine
             stella_config_name = "Helias 4"
             ! Reference point where all the other variables are determined from
             ! Plasma outer radius
             stella_config_rmajor_ref = 17.6D0
             stella_config_rminor_ref = 2.0D0
             stella_config_aspect_ref =  8.8D0

             ! Coil radii
             stella_config_coil_rmajor = 18.39D0
             stella_config_coil_rminor = 4.94D0

             stella_config_bt_ref = 5.6D0
             stella_config_WP_area = 0.8d0*0.6d0
             stella_config_WP_bmax = 11.51d0

             stella_config_symmetry = 4
             stella_config_coilspermodule = 10
             stella_config_a1 = 0.676D0
             stella_config_a2 = 0.029D0
             stella_config_plasma_volume =   1380.0D0
             stella_config_dmin = 1.08D0
             stella_config_max_portsize_width = 3.24D0

             stella_config_plasma_surface = 1900.0D0
             stella_config_maximal_coil_height = 13.34D0  ! [m] Full height max point to min point

             stella_config_coilsurface =  4100.0D0! Coil surface, dimensionfull. At reference point

             stella_config_coillength = 1435.07D0 ! Central filament length of machine with outer radius 1m.

             stella_config_I0 = 13.146D0 ! Coil Current needed to produce b0 on axis in [MA] at reference point
             stella_config_inductance = 1290.4D-6 ! inductance/R*A^2 in muH

             stella_config_WP_ratio = 1.3D0

             stella_config_max_force_density = 120.0d0 ! [MN/m^3]
             stella_config_max_force_density_MNm = 98.0d0 ! [MN/m]


             stella_config_max_lateral_force_density = 87.9d0 ! [MN/m^3]
             stella_config_max_radial_force_density = 109.9d0   ! [MN/m^3]

             stella_config_centering_force_max_MN = 226.0d0
             stella_config_centering_force_min_MN = -35.3d0
             stella_config_centering_force_avg_MN = 125.8d0

             stella_config_min_plasma_coil_distance = 1.7d0
             stella_config_derivative_min_LCFS_coils_dist = -1.0d0 ! this is approximated for now

             stella_config_min_bend_radius = 0.86d0 ! [m]

             stella_config_neutron_peakfactor = 1.6d0

             stella_config_epseff = 0.015d0


          case(3)
             ! Helias 3 Machine
             stella_config_name = "Helias 3"
             ! Reference point where all the other variables are determined from
             ! Plasma outer radius
             stella_config_rmajor_ref = 13.86d0
             stella_config_rminor_ref = 2.18d0
             stella_config_aspect_ref =  6.36d0

             ! Coil radii
             stella_config_coil_rmajor = 14.53D0
             stella_config_coil_rminor = 6.12D0

             stella_config_bt_ref = 5.6D0
             stella_config_WP_bmax = 12.346d0
             stella_config_WP_area = 0.8d0*0.6d0

             stella_config_symmetry = 3
             stella_config_coilspermodule = 10

             ! Bmax fit parameters
             stella_config_a1 = 0.56D0
             stella_config_a2 = 0.030D0

             stella_config_plasma_volume =   1300.8D0
             stella_config_dmin = 1.145D0
             stella_config_max_portsize_width = 3.24D0 !??? guess. not ready yet

             stella_config_plasma_surface = 1600.00D0

             stella_config_maximal_coil_height = 17.74D0! [m] Full height max point to min point

             stella_config_coilsurface = 4240.0D0 ! Coil surface, dimensionfull. At reference point

             stella_config_coillength = 1287.3D0 ! Central filament length of machine with outer radius 1m.

             stella_config_I0 = 14.23D0 ! Coil Current needed to produce 1T on axis in [MA] at outer radius 1m
             stella_config_inductance = 1250.7D-6 ! inductance in muH

             stella_config_WP_ratio = 1.3D0

             stella_config_max_force_density = 120.0d0
             stella_config_max_force_density_MNm = 98.0d0 ! [MN/m]


             stella_config_max_lateral_force_density = 96.6d0 ! [MN/m^3]
             stella_config_max_radial_force_density = 130.5d0   ! [MN/m^3]

             stella_config_centering_force_max_MN = 428.1d0
             stella_config_centering_force_min_MN = -70.3d0
             stella_config_centering_force_avg_MN = 240.9d0

             stella_config_min_plasma_coil_distance = 1.78d0
             stella_config_derivative_min_LCFS_coils_dist = -1.0d0 ! this is approximated for now

             stella_config_min_bend_radius = 1.145d0 ! [m]

             stella_config_neutron_peakfactor = 1.6d0

             stella_config_epseff = 0.015d0


          case(4)
             ! w7x30 Machine
             stella_config_name = "W7X-30"
             ! Reference point where all the other variables are determined from
             ! Plasma outer radius
             stella_config_rmajor_ref = 5.50D0
             stella_config_rminor_ref = 0.49D0
             stella_config_aspect_ref =  11.2D0

             ! Coil radii
             stella_config_coil_rmajor = 5.62D0
             stella_config_coil_rminor = 1.36D0


             stella_config_bt_ref = 3.0D0
             stella_config_WP_area = 0.18d0*0.15d0
             stella_config_WP_bmax = 10.6d0

             stella_config_symmetry = 5
             stella_config_coilspermodule = 6
             stella_config_a1 = 0.98D0
             stella_config_a2 = 0.041D0
             stella_config_plasma_volume =   26.4D0
             stella_config_dmin = 0.21D0
             stella_config_max_portsize_width = 0.5D0

             stella_config_plasma_surface = 128.3D0
             stella_config_maximal_coil_height = 3.6D0  ! [m] Full height max point to min point

             stella_config_coilsurface =  370.0D0! Coil surface, dimensionfull. At reference point

             stella_config_coillength = 303.4D0 ! Central filament length of machine with outer radius 1m.

             stella_config_I0 = 2.9D0 ! Coil Current needed to produce b0 on axis in [MA] at reference point
             stella_config_inductance = 252.7D-6 ! inductance/R*A^2 in muH

             stella_config_WP_ratio = 1.2D0

             stella_config_max_force_density = 350.0d0 ! [MN/m^3]
             stella_config_max_force_density_MNm = 98.0d0 ! [MN/m]


             stella_config_max_lateral_force_density = 271.1d0 ! [MN/m^3]
             stella_config_max_radial_force_density = 305.2d0   ! [MN/m^3]

             stella_config_centering_force_max_MN = 7.95d0
             stella_config_centering_force_min_MN = -2.15d0
             stella_config_centering_force_avg_MN = 3.46d0

             stella_config_min_plasma_coil_distance = 0.45D0
             stella_config_derivative_min_LCFS_coils_dist = -1.0d0 ! this is approximated for now

             stella_config_min_bend_radius = 0.186d0 ! [m]

             stella_config_neutron_peakfactor = 1.6d0

             stella_config_epseff = 0.015d0

          case(5)
             ! w7x50 Machine
             stella_config_name = "W7X-50"
             ! Reference point where all the other variables are determined from
             ! Plasma outer radius
             stella_config_rmajor_ref = 5.5D0
             stella_config_rminor_ref = 0.49D0
             stella_config_aspect_ref =   11.2D0

             ! Coil radii
             stella_config_coil_rmajor = 5.62d0
             stella_config_coil_rminor = 1.18D0


             stella_config_bt_ref = 3.0D0
             stella_config_WP_area = 0.18d0*0.15d0
             stella_config_WP_bmax = 6.3d0

             stella_config_symmetry = 5
             stella_config_coilspermodule = 10
             stella_config_a1 = 0.66D0
             stella_config_a2 = 0.025D0
             stella_config_plasma_volume =   26.4D0
             stella_config_dmin = 0.28D0
             stella_config_max_portsize_width = 0.3D0

             stella_config_plasma_surface = 128.3D0
             stella_config_maximal_coil_height = 3.1D0  ! [m] Full height max point to min point

             stella_config_coilsurface =  299.85D0! Coil surface, dimensionfull. At reference point

             stella_config_coillength = 420.67D0 ! Central filament length of machine with outer radius 1m.

             stella_config_I0 = 1.745D0 ! Coil Current needed to produce b0 on axis in [MA] at reference point
             stella_config_inductance = 412.4D-6 ! inductance/R*A^2 in muH

             stella_config_WP_ratio = 1.2D0

             stella_config_max_force_density = 250.0d0 ! [MN/m^3]
             stella_config_max_force_density_MNm = 98.0d0 ! [MN/m]


             stella_config_max_lateral_force_density = 116.4d0 ! [MN/m^3]
             stella_config_max_radial_force_density = 148.d0   ! [MN/m^3]

             stella_config_centering_force_max_MN = 2.99d0
             stella_config_centering_force_min_MN = -1.29d0
             stella_config_centering_force_avg_MN = 1.61d0

             stella_config_min_plasma_coil_distance = 0.39D0
             stella_config_derivative_min_LCFS_coils_dist = -1.0d0 ! this is approximated for now

             stella_config_min_bend_radius = 0.39d0 ! [m]

             stella_config_neutron_peakfactor = 1.6d0

             stella_config_epseff = 0.015d0



          case(6)
             ! Init from json
             ! This requires a file called stella_config.json in the working directory.
             ! It can either be prepared manually or it can be produced automatically based on a VMEC netcdf
             ! file and a coils file
             ! by the pre-sPROCESS Code, https://gitlab.mpcdf.mpg.de/jtl/sprocess/ by jorrit.lion@ipp.mpg.de
             call stella_config_json()
          case default
             ! Return some error here. The index is not implemented yet.
             write(*,*)'ERROR in initialization of stellarator config. No such istell: ',index
       end select

    end subroutine new_stella_config



    subroutine stella_config_json()

        !! Initialises the effective stellarator values using a json input file
        !! author: J Lion, IPP Greifswald
        !! None
        !! This routine reads in all effective variables that
        !! are given needed by the 'constructor' stella_config
        !! <P>The effective values are read in from a JSON-format file.
        !! None
        !
        ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        use fson_library, only: fson_parse, fson_value, fson_get, fson_destroy
        use global_variables, only: output_prefix



        !  Arguments

        !  Local variables

        integer :: n_values
        character(len=180) :: filename
        type(fson_value), pointer :: stellafile

        real(dp), dimension(:), allocatable :: nustar,d11,d13
        !type(stella_config), allocatable, dimension(:) :: stella_json

        ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        !  Parse the json file

        filename = trim(output_prefix) // 'stella_conf.json'
        stellafile => fson_parse(trim(filename))

        !  Extract information arrays from the file

        call fson_get(stellafile, "name", stella_config_name)
        call fson_get(stellafile, "symmetry", stella_config_symmetry)

        call fson_get(stellafile, "coilspermodule", stella_config_coilspermodule)
        call fson_get(stellafile, "rmajor_ref", stella_config_rmajor_ref)
        call fson_get(stellafile, "rminor_ref", stella_config_rminor_ref)
        call fson_get(stellafile, "coil_rmajor", stella_config_coil_rmajor)
        call fson_get(stellafile, "coil_rminor", stella_config_coil_rminor)
        call fson_get(stellafile, "aspect_ref", stella_config_aspect_ref)
        call fson_get(stellafile, "bt_ref", stella_config_bt_ref)
        call fson_get(stellafile, "WP_area", stella_config_WP_area)
        call fson_get(stellafile, "WP_bmax", stella_config_WP_bmax)
        call fson_get(stellafile, "i0", stella_config_i0)
        call fson_get(stellafile, "a1", stella_config_a1)
        call fson_get(stellafile, "a2", stella_config_a2)
        call fson_get(stellafile, "dmin", stella_config_dmin)
        call fson_get(stellafile, "inductance", stella_config_inductance)
        call fson_get(stellafile, "coilsurface", stella_config_coilsurface)
        call fson_get(stellafile, "coillength", stella_config_coillength)
        call fson_get(stellafile, "max_portsize_width", stella_config_max_portsize_width)
        call fson_get(stellafile, "maximal_coil_height", stella_config_maximal_coil_height)
        call fson_get(stellafile, "min_plasma_coil_distance", stella_config_min_plasma_coil_distance)
        call fson_get(stellafile, "derivative_min_LCFS_coils_dist", stella_config_derivative_min_LCFS_coils_dist)
        call fson_get(stellafile, "plasma_volume", stella_config_plasma_volume)
        call fson_get(stellafile, "plasma_surface", stella_config_plasma_surface)
        call fson_get(stellafile, "WP_ratio", stella_config_WP_ratio)
        call fson_get(stellafile, "max_force_density", stella_config_max_force_density)
        call fson_get(stellafile, "max_force_density_MNm", stella_config_max_force_density_MNm)
        call fson_get(stellafile, "min_bend_radius", stella_config_min_bend_radius)
        call fson_get(stellafile, "epseff", stella_config_epseff)

        call fson_get(stellafile, "max_lateral_force_density", stella_config_max_lateral_force_density)
        call fson_get(stellafile, "max_radial_force_density", stella_config_max_radial_force_density)

        call fson_get(stellafile, "centering_force_max_MN", stella_config_centering_force_max_MN)
        call fson_get(stellafile, "centering_force_min_MN", stella_config_centering_force_min_MN)
        call fson_get(stellafile, "centering_force_avg_MN", stella_config_centering_force_avg_MN)

        call fson_get(stellafile, "neutron_peakfactor", stella_config_neutron_peakfactor)


        call fson_get(stellafile, "number_nu_star", n_values)

        if (allocated(sc_D11_star_mono_input)) deallocate(sc_D11_star_mono_input)
        if (allocated(sc_nu_star_mono_input)) deallocate(sc_nu_star_mono_input)
        allocate(sc_D11_star_mono_input(n_values))
        allocate(sc_nu_star_mono_input(n_values))


        call fson_get(stellafile, "D11_star_mono_input", sc_D11_star_mono_input)
        call fson_get(stellafile, "nu_star_mono_input", sc_nu_star_mono_input)




        !  Clean up
        call fson_destroy(stellafile)

    end subroutine stella_config_json


    subroutine stella_error(index,keyname)
    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Gives the errors of a the stellarator module

    integer, intent(in) :: index
    character, intent(in) :: keyname


    select case (index)
    case(1)
        ! Error reading in a json attribute
        ! Not used yet because I don't know how.
        write(*,*)'ERROR in initialization of stellarator config. Missing json key: ',keyname


    case default
        ! Return some error here. The index is not implemented yet.
        write(*,*)'ERROR in stella_error! No such error index: ',index

    end select


    end subroutine stella_error


 end module stellarator_configuration

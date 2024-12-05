module neoclassics_constants
    integer, parameter :: no_roots = 30 ! Number of Gauss laguerre roots
end module neoclassics_constants

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module neoclassics_module

    !! Module containing neoclassical computations
    !! author: J Lion, IPP Greifswald
    !! Formulas used are described in:
    !! Beidler (2013), https://doi.org/10.1088/0029-5515/51/7/076001
    !
    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#ifndef dp
    use, intrinsic :: iso_fortran_env, only: dp=>real64
#endif

    use neoclassics_constants, only: no_roots

    implicit none

    public



    ! These module variables were taken from the previous Neoclassics derived type. This derived type contained
    ! nested derived types 'Gauss_leguerre' and 'Profile_values'. They were made module variables for ease
    ! of import into the stellarator module such that Python conversion could take place.
    ! S Gubbins 12/09/2022

    character, dimension(4) :: species = (/"e","D","T","a"/)
    !  Species that are considered
    real(dp), dimension(4) :: densities
    !  Densities of the species that are considered [/m3]
    real(dp), dimension(4) :: temperatures
    !  Temperature of the species that are considered [J]
    real(dp), dimension(4) :: dr_densities
    !  Radial derivative of the density of the species [/m3]
    real(dp), dimension(4) :: dr_temperatures
    !  Radial derivative of the temperature of the species [J]
    real(dp), dimension(no_roots) :: roots = 0
    !  Gauss Laguerre Roots
    real(dp), dimension(no_roots) :: weights = 0
    !  Gauss Laguerre Weights
    real(dp), dimension(4,no_roots) :: nu = 0
    !  90-degree deflection frequency on GL roots
    real(dp), dimension(4,no_roots) :: nu_star = 0
    !  Dimensionless deflection frequency
    real(dp), dimension(4) :: nu_star_averaged = 0
    !  Maxwellian averaged dimensionless 90-degree deflection frequency for electrons (index 1) and ions (index 2)
    real(dp), dimension(4,no_roots) :: vd = 0
    !  Drift velocity on GL roots
    real(dp), dimension(4,no_roots) :: KT = 0
    !  Thermal energy on GL roots
    real(dp) :: Er = 0.0
    !  Radial electrical field [V/m]
    real(dp) :: iota = 1.0d0
    !  Iota (1/safety factor)
    real(dp), dimension(4,no_roots) :: D11_mono = 0
    !  Radial monoenergetic transport coefficient on GL roots (species dependent)
    real(dp), dimension(4,no_roots) :: D11_plateau = 0
    !  Toroidal monoenergetic transport coefficient as given by the stellarator
    !  input json file as function of nu_star, normalized by the banana value.
    real(dp), dimension(4) :: D111 = 0
    !  Radial integrated transport coefficient (n=1) (species dependent)
    real(dp), dimension(4) :: D112 = 0
    !  Radial integrated transport coefficient (n=2) (species dependent)
    real(dp), dimension(4) :: D113 = 0
    !  Radial integrated transport coefficient (n=3) (species dependent)
    real(dp), dimension(4) :: q_flux = 0
    !  energy transport flux (J/m2)
    real(dp), dimension(4) :: Gamma_flux = 0
    !  energy flux from particle transport
    real(dp), dimension(no_roots) :: D31_mono = 0
    !  Toroidal monoenergetic transport coefficient
    real(dp) :: eps_eff = 1d-5
    !  Epsilon effective (used in neoclassics_calc_D11_mono)
    real(dp) :: r_eff = 0




contains
    subroutine init_neoclassics_module
        !! Initialise module variables
        implicit none
        species = (/"e","D","T","a"/)
        densities = 0
        temperatures = 0
        dr_densities = 0
        dr_temperatures =0
        roots = 0
        weights = 0
        nu = 0
        nu_star = 0
        nu_star_averaged = 0
        vd = 0
        KT = 0
        Er = 0.0
        iota = 1.0d0
        D11_mono = 0
        D11_plateau = 0
        D111 = 0
        D112 = 0
        D113 = 0
        q_flux = 0
        Gamma_flux = 0
        D31_mono = 0
        eps_eff = 1d-5
    end subroutine init_neoclassics_module

end module neoclassics_module

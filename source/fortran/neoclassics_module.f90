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
    !  Radial monoenergetic transport coefficient on GL roots (species dependent)
    real(dp), dimension(:), allocatable :: nu_star_mono_input
    !  Radial monoenergetic transport coefficient as given by the stellarator input json
    !  on GL roots (species dependent)
    real(dp), dimension(:), allocatable :: D11_star_mono_input
    !  Radial monoenergetic transport coefficient as given by the stellarator input json
    !  as function of nu_star, normalized by the plateau value.
    real(dp), dimension(:), allocatable :: D13_star_mono_input
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
        nu_star_mono_input = 0
        D11_star_mono_input = 0
        D13_star_mono_input = 0
        D111 = 0
        D112 = 0
        D113 = 0
        q_flux = 0
        Gamma_flux = 0
        D31_mono = 0
        eps_eff = 1d-5
    end subroutine init_neoclassics_module

    subroutine init_neoclassics(r_eff, eps_eff, iota)
    !! Constructor of the neoclassics object from the effective radius,
    !! epsilon effective and iota only.
    !
    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    real(dp), intent(in) :: r_eff
    real(dp), dimension(4,no_roots) :: mynu
    real(dp), intent(in)::eps_eff, iota

    !! This should be called as the standard constructor
    call init_profile_values_from_PROCESS(r_eff, densities, temperatures, dr_densities, dr_temperatures)
    call gauss_laguerre_30_roots(roots)
    call gauss_laguerre_30_weights(weights)


    KT = neoclassics_calc_KT()
    nu = neoclassics_calc_nu()
    nu_star = neoclassics_calc_nu_star()
    nu_star_averaged = neoclassics_calc_nu_star_fromT(iota)
    vd = neoclassics_calc_vd()

    D11_plateau = neoclassics_calc_D11_plateau()

    D11_mono = neoclassics_calc_D11_mono(eps_eff) !for using epseff

    !alternatively use:  = myneo%interpolate_D11_mono() !

    D111 = neoclassics_calc_D111()

    D112 = neoclassics_calc_D112()
    D113 = neoclassics_calc_D113()

    Gamma_flux = neoclassics_calc_Gamma_flux(densities, temperatures, dr_densities, dr_temperatures)
    q_flux = neoclassics_calc_q_flux()

    ! Return:

    end subroutine init_neoclassics


    function neoclassics_calc_KT() result(KK)
        !! Calculates the energy on the given grid
        !! which is given by the gauss laguerre roots.
        !
        ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        use const_and_precisions, only: e_,c_,me_,mp_,keV_

        real(dp), dimension(no_roots) ::  K
        real(dp), dimension(4,no_roots) :: KK

        K = roots/keV_

        KK(1,:) = K * temperatures(1) ! electrons
        KK(2,:) = K * temperatures(2) ! deuterium
        KK(3,:) = K * temperatures(3) ! tritium
        KK(4,:) = K * temperatures(4) ! helium

    end function neoclassics_calc_KT

    function neoclassics_calc_Gamma_flux(densities, temperatures, dr_densities, dr_temperatures)
        !! Calculates the Energy flux by neoclassical particle transport
        !
        ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        real(dp),dimension(4) :: neoclassics_calc_Gamma_flux, densities, dr_densities, z, temperatures, dr_temperatures


        z = (/-1.0,1.0,1.0,2.0/)

        neoclassics_calc_Gamma_flux = - densities * D111 * ((dr_densities/densities - z * Er/temperatures)+ &
                        (D112/D111-3.0/2.0) * dr_temperatures/temperatures )

    end function neoclassics_calc_Gamma_flux

    function neoclassics_calc_q_flux()
        !! Calculates the Energy flux by neoclassicsal energy transport
        !
        ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        real(dp),dimension(4) ::  z, neoclassics_calc_q_flux


        ! densities = densities
        ! temps = temperatures
        ! dr_densities = dr_densities
        ! dr_temps = dr_temperatures

        z = (/-1.0,1.0,1.0,2.0/)

        q_flux = - densities * temperatures * D112 * ((dr_densities/densities - z * Er/temperatures) + &
                        (D113/D112-3.0/2.0) * dr_temperatures/temperatures )

        neoclassics_calc_q_flux = q_flux
    end function neoclassics_calc_q_flux

    function neoclassics_calc_D11_mono(eps_eff) result(D11_mono)
        !! Calculates the monoenergetic radial transport coefficients
        !! using epsilon effective.
        !
        ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        use const_and_precisions, only: pi

        real(dp),dimension(4,no_roots) :: D11_mono
        real(dp), intent(in):: eps_eff

        D11_mono = 4.0d0/(9.0d0*pi) * (2.0d0 * eps_eff)**(3.0d0/2.0d0) &
                    * vd**2/nu

    end function neoclassics_calc_D11_mono

    function neoclassics_calc_D11_plateau() result(D11_plateau)
        !! Calculates the plateau transport coefficients (D11_star sometimes)
        !
        ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        use const_and_precisions, only: pi, me_, mp_, c_
        use physics_variables, only: rmajor

        real(dp),dimension(4,no_roots) :: D11_plateau, v
        real(dp),dimension(4) :: mass

        mass = (/me_,mp_*2.0d0,mp_*3.0d0,mp_*4.0d0/)

        v(1,:) = c_ * sqrt(1.0d0-(KT(1,:)/(mass(1) * c_**2)+1)**(-1))
        v(2,:) = c_ * sqrt(1.0d0-(KT(2,:)/(mass(2) * c_**2)+1)**(-1))
        v(3,:) = c_ * sqrt(1.0d0-(KT(3,:)/(mass(3) * c_**2)+1)**(-1))
        v(4,:) = c_ * sqrt(1.0d0-(KT(4,:)/(mass(4) * c_**2)+1)**(-1))

        D11_plateau = pi/4.0 * vd**2 * rmajor/ iota / v

    end function neoclassics_calc_D11_plateau

    function neoclassics_interpolate_D11_mono() result(D11_mono)
        !! Interpolates the D11 coefficients on the Gauss laguerre grid
        !! (This method is unused as of now, but is needed when taking D11 explicitely as input)
        !
        ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        use const_and_precisions, only: pi
        use maths_library, only: find_y_nonuniform_x
        ! use grad_func, only: interp1_ef

        real(dp),dimension(4,no_roots) :: D11_mono
        integer :: ii,jj

        do ii = 1,4
            do jj = 1,no_roots
                D11_mono(ii,jj) = find_y_nonuniform_x(nu_star(ii,jj),nu_star_mono_input, &
                                                      D11_star_mono_input,size(nu_star_mono_input)) * &
                                  D11_plateau(ii,jj)
            end do
        end do

    end function neoclassics_interpolate_D11_mono

    function neoclassics_calc_vd()
        !! Calculates the drift velocities
        !
        ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        use const_and_precisions, only: e_
        use physics_variables, only: rmajor, bt

        real(dp), dimension(no_roots) :: vde,vdT,vdD,vda, K
        real(dp), dimension(4,no_roots) :: vd,neoclassics_calc_vd

        K = roots

        vde = K * temperatures(1)/(e_ * rmajor * bt)
        vdD = K * temperatures(2)/(e_ * rmajor * bt)
        vdT = K * temperatures(3)/(e_ * rmajor * bt)
        vda = K * temperatures(4)/(2.0*e_ * rmajor * bt)

        vd(1,:) = vde
        vd(2,:) = vdD
        vd(3,:) = vdT
        vd(4,:) = vda

        neoclassics_calc_vd = vd
    end function neoclassics_calc_vd

    function neoclassics_calc_nu_star() result(nu_star)
        !! Calculates the normalized collision frequency
        !
        ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        use const_and_precisions, only: e_,c_,me_,mp_
        use physics_variables, only: rmajor

        real(dp), dimension(no_roots) ::  K
        real(dp), dimension(4,no_roots) :: v,nu_star,KK
        real(dp), dimension(4) :: mass

        K = roots

        KK(1,:) = K * temperatures(1)
        KK(2,:) = K * temperatures(2)
        KK(3,:) = K * temperatures(3)
        KK(4,:) = K * temperatures(4)

        mass = (/me_,mp_*2.0d0,mp_*3.0d0,mp_*4.0d0/)

        v(1,:) = c_ * sqrt(1.0d0-(KK(1,:)/(mass(1) * c_**2)+1)**(-1))
        v(2,:) = c_ * sqrt(1.0d0-(KK(2,:)/(mass(2) * c_**2)+1)**(-1))
        v(3,:) = c_ * sqrt(1.0d0-(KK(3,:)/(mass(3) * c_**2)+1)**(-1))
        v(4,:) = c_ * sqrt(1.0d0-(KK(4,:)/(mass(4) * c_**2)+1)**(-1))

        nu_star = rmajor * nu/(iota*v)

    end function neoclassics_calc_nu_star


    function neoclassics_calc_nu_star_fromT(iota)
        !! Calculates the collision frequency
        !
        ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        use const_and_precisions, only: pi, me_, mp_, eps0_,e_, keV_
        use physics_variables, only: rmajor, te,ti, dene,deni, dnalp, fdeut

        real(dp),dimension(4) :: neoclassics_calc_nu_star_fromT
        real(dp) :: t,erfn,phixmgx,expxk,xk, lnlambda,x,v
        real(dp),dimension(4) :: temp, mass,density,z
        real(dp) :: iota

        integer :: jj,kk

        temp = (/te,ti,ti,ti /) * keV_
        density = (/dene,deni * fdeut,deni*(1-fdeut),dnalp /)

        !          e      D      T         a (He)
        mass = (/me_,mp_*2.0d0,mp_*3.0d0,mp_*4.0d0/)
        z = (/-1.0d0,1.0d0,1.0d0,2.0d0/) * e_

        ! transform the temperature back in eV
        ! Formula from L. Spitzer.Physics of fully ionized gases.  Interscience, New York, 1962
        lnlambda = 32.2d0 - 1.15d0*log10(density(1)) + 2.3d0*log10(temp(1)/e_)

        neoclassics_calc_nu_star_fromT(:) = 0.0d0

        do jj = 1, 4
            v = sqrt(2d0 * temp(jj)/mass(jj))
            do kk = 1,4
                xk = (mass(kk)/mass(jj))*(temp(jj)/temp(kk))

                if (xk < 200.d0) then
                    expxk = exp(-xk)
                else
                    expxk = 0.0d0
                endif

                t = 1.0d0/(1.0d0+0.3275911d0*sqrt(xk))
                erfn = 1.0d0-t*(.254829592d0 + t*(-.284496736d0 + t*(1.421413741d0       &
                        + t*(-1.453152027d0 +t*1.061405429d0))))*expxk
                phixmgx = (1.0d0-0.5d0/xk)*erfn + expxk/sqrt(pi*xk)
                neoclassics_calc_nu_star_fromT(jj) = neoclassics_calc_nu_star_fromT(jj) + density(kk)*(z(jj)*z(kk))**2 &
                            *lnlambda*phixmgx/(4.0d0*pi*eps0_**2*mass(jj)**2*v**4) * rmajor/iota
            enddo
        enddo

    end function neoclassics_calc_nu_star_fromT

    function neoclassics_calc_D111()
        !! Calculates the integrated radial transport coefficients (index 1)
        !! It uses Gauss laguerre integration
        !! https://en.wikipedia.org/wiki/Gauss%E2%80%93Laguerre_quadrature
        !
        ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        use const_and_precisions, only: pi

        real(dp),dimension(4) :: D111, neoclassics_calc_D111

        real(dp),dimension(no_roots) :: xi,wi

        xi = roots
        wi = weights

        D111(1) = sum(2.0d0/sqrt(pi) * D11_mono(1,:) * xi**(1.0d0-0.5d0) * wi)
        D111(2) = sum(2.0d0/sqrt(pi) * D11_mono(2,:) * xi**(1.0d0-0.5d0) * wi)
        D111(3) = sum(2.0d0/sqrt(pi) * D11_mono(3,:) * xi**(1.0d0-0.5d0) * wi)
        D111(4) = sum(2.0d0/sqrt(pi) * D11_mono(4,:) * xi**(1.0d0-0.5d0) * wi)

        neoclassics_calc_D111 = D111

    end function neoclassics_calc_D111

    function neoclassics_calc_D112()
        !! Calculates the integrated radial transport coefficients (index 2)
        !! It uses Gauss laguerre integration
        !! https://en.wikipedia.org/wiki/Gauss%E2%80%93Laguerre_quadrature
        !
        ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        use const_and_precisions, only: pi

        real(dp),dimension(4) :: D112, neoclassics_calc_D112

        real(dp),dimension(no_roots) :: xi,wi

        xi = roots
        wi = weights

        D112(1) = sum(2.0d0/sqrt(pi) * D11_mono(1,:) * xi**(2.0d0-0.5d0) * wi)
        D112(2) = sum(2.0d0/sqrt(pi) * D11_mono(2,:) * xi**(2.0d0-0.5d0) * wi)
        D112(3) = sum(2.0d0/sqrt(pi) * D11_mono(3,:) * xi**(2.0d0-0.5d0) * wi)
        D112(4) = sum(2.0d0/sqrt(pi) * D11_mono(4,:) * xi**(2.0d0-0.5d0) * wi)

        neoclassics_calc_D112 = D112

    end function neoclassics_calc_D112

    function neoclassics_calc_D113()
        !! Calculates the integrated radial transport coefficients (index 3)
        !! It uses Gauss laguerre integration
        !! https://en.wikipedia.org/wiki/Gauss%E2%80%93Laguerre_quadrature
        !
        ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        use const_and_precisions, only: pi

        real(dp),dimension(4) :: D113, neoclassics_calc_D113

        real(dp),dimension(no_roots) :: xi,wi

        xi = roots
        wi = weights

        D113(1) = sum(2.0d0/sqrt(pi) * D11_mono(1,:) * xi**(3.0d0-0.5d0) * wi)
        D113(2) = sum(2.0d0/sqrt(pi) * D11_mono(2,:) * xi**(3.0d0-0.5d0) * wi)
        D113(3) = sum(2.0d0/sqrt(pi) * D11_mono(3,:) * xi**(3.0d0-0.5d0) * wi)
        D113(4) = sum(2.0d0/sqrt(pi) * D11_mono(4,:) * xi**(3.0d0-0.5d0) * wi)

        neoclassics_calc_D113 = D113
    end function neoclassics_calc_D113

    function neoclassics_calc_nu()
        !! Calculates the collision frequency
        !
        ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        use const_and_precisions, only: pi, me_, mp_, eps0_,e_

        real(dp),dimension(4,no_roots) :: neoclassics_calc_nu
        real(dp) :: t,erfn,phixmgx,expxk,xk, lnlambda,x,v
        real(dp),dimension(4) :: temp, mass,density,z

        integer :: jj,ii,kk

        temp = temperatures
        density = densities

        !          e      D      T         a (He)
        mass = (/me_,mp_*2.0d0,mp_*3.0d0,mp_*4.0d0/)
        z = (/-1.0d0,1.0d0,1.0d0,2.0d0/) * e_

        ! transform the temperature back in eV
        ! Formula from L. Spitzer.Physics of fully ionized gases.  Interscience, New York, 1962
        lnlambda = 32.2d0 - 1.15d0*log10(density(1)) + 2.3d0*log10(temp(1)/e_)

        neoclassics_calc_nu(:,:) = 0.0

        do jj = 1, 4
           do ii = 1, no_roots
              x = roots(ii)
              do kk = 1,4
                 xk = (mass(kk)/mass(jj))*(temp(jj)/temp(kk))*x
                 expxk = exp(-xk)
                 t = 1.0d0/(1.0d0+0.3275911d0*sqrt(xk))
                 erfn = 1.0d0-t*(.254829592d0 + t*(-.284496736d0 + t*(1.421413741d0       &
                         + t*(-1.453152027d0 +t*1.061405429d0))))*expxk
                 phixmgx = (1.0d0-0.5d0/xk)*erfn + expxk/sqrt(pi*xk)
                 v = sqrt(2.0d0*x*temp(jj)/mass(jj))
                 neoclassics_calc_nu(jj,ii) = neoclassics_calc_nu(jj,ii) + density(kk)*(z(jj)*z(kk))**2 &
                              *lnlambda *phixmgx/(4.0*pi*eps0_**2*mass(jj)**2*v**3)
              enddo
           enddo
        enddo

    end function neoclassics_calc_nu

    subroutine init_profile_values_from_PROCESS(rho, densities, temperatures, dr_densities, dr_temperatures)
        !! Initializes the profile_values object from PROCESS' parabolic profiles
        !
        ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        use physics_variables, only: ne0,te0,alphan,&
                                     alphat,ti0,ni0,fdeut, dnalp, rminor
        use const_and_precisions, only: keV_

        real(dp), intent(in) :: rho

        real(dp),dimension(4) :: dens,temp, dr_dens, dr_temp
        real(dp) :: dense, densD,densT,densa, &
                    tempD,tempT,tempa,tempe, &
                    dr_tempe, dr_tempT, dr_tempD, dr_tempa,&
                    dr_dense, dr_densT, dr_densD, dr_densa, r
        real(dp), dimension(4), intent(out):: densities
        real(dp), dimension(4), intent(out):: temperatures
        real(dp), dimension(4), intent(out):: dr_densities
        real(dp), dimension(4), intent(out):: dr_temperatures

        r = rho * rminor

        tempe = te0 * (1-rho**2)**alphat * keV_ ! To SI units bc.. convenience I guess?
        tempT = ti0 * (1-rho**2)**alphat * keV_
        tempD = ti0 * (1-rho**2)**alphat * keV_
        tempa = ti0 * (1-rho**2)**alphat * keV_

        dense = ne0 * (1-rho**2)**alphan
        densT = (1-fdeut) * ni0 * (1-rho**2)**alphan
        densD = fdeut *ni0 * (1-rho**2)**alphan
        densa = dnalp*(1+alphan) * (1-rho**2)**alphan

        ! Derivatives in real space
        dr_tempe = -2.0d0 * 1.0d0/rminor * te0 * rho * (1.0d0-rho**2)**(alphat-1.0d0) * alphat * keV_
        dr_tempT = -2.0d0 * 1.0d0/rminor * ti0 * rho * (1.0d0-rho**2)**(alphat-1.0d0) * alphat * keV_
        dr_tempD = -2.0d0 * 1.0d0/rminor * ti0 * rho * (1.0d0-rho**2)**(alphat-1.0d0) * alphat * keV_
        dr_tempa = -2.0d0 * 1.0d0/rminor * ti0 * rho * (1.0d0-rho**2)**(alphat-1.0d0) * alphat * keV_

        dr_dense = -2.0d0 * 1.0d0/rminor * rho * ne0 *             (1.0d0-rho**2)**(alphan-1.0d0) * alphan
        dr_densT = -2.0d0 * 1.0d0/rminor * rho * (1-fdeut) * ni0 * (1.0d0-rho**2)**(alphan-1.0d0) * alphan
        dr_densD = -2.0d0 * 1.0d0/rminor * rho * fdeut *ni0 *      (1.0d0-rho**2)**(alphan-1.0d0) * alphan
        dr_densa = -2.0d0 * 1.0d0/rminor * rho * dnalp*(1+alphan)* (1.0d0-rho**2)**(alphan-1.0d0) * alphan

        dens(1) = dense
        dens(2) = densD
        dens(3) = densT
        dens(4) = densa

        temp(1) = tempe
        temp(2) = tempD
        temp(3) = tempT
        temp(4) = tempa

        dr_dens(1) = dr_dense
        dr_dens(2) = dr_densD
        dr_dens(3) = dr_densT
        dr_dens(4) = dr_densa

        dr_temp(1) = dr_tempe
        dr_temp(2) = dr_tempD
        dr_temp(3) = dr_tempT
        dr_temp(4) = dr_tempa

        densities = dens
        temperatures = temp
        dr_densities = dr_dens
        dr_temperatures = dr_temp

    end subroutine init_profile_values_from_PROCESS

    subroutine gauss_laguerre_30_roots(roots)

        !! Sets the gauss Laguerre roots and weights for 30
        !! discretization points. Used for integration in this module.
        !! roots
        !
        ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        real(dp), dimension(no_roots), intent(out):: roots
        roots= (/4.740718054080526184d-02,&
                                    2.499239167531593919d-01,&
                                    6.148334543927683749d-01,&
                                    1.143195825666101451d+00,&
                                    1.836454554622572344d+00,&
                                    2.696521874557216147d+00,&
                                    3.725814507779509288d+00,&
                                    4.927293765849881879d+00,&
                                    6.304515590965073635d+00,&
                                    7.861693293370260349d+00,&
                                    9.603775985479263255d+00,&
                                    1.153654659795613924d+01,&
                                    1.366674469306423489d+01,&
                                    1.600222118898106771d+01,&
                                    1.855213484014315029d+01,&
                                    2.132720432178312819d+01,&
                                    2.434003576453269346d+01,&
                                    2.760555479678096091d+01,&
                                    3.114158670111123683d+01,&
                                    3.496965200824907072d+01,&
                                    3.911608494906788991d+01,&
                                    4.361365290848483056d+01,&
                                    4.850398616380419980d+01,&
                                    5.384138540650750571d+01,&
                                    5.969912185923549686d+01,&
                                    6.618061779443848991d+01,&
                                    7.344123859555988076d+01,&
                                    8.173681050672767867d+01,&
                                    9.155646652253683726d+01,&
                                    1.041575244310588886d+02/)

    end subroutine gauss_laguerre_30_roots

    subroutine gauss_laguerre_30_weights(weights)

        real(dp), dimension(no_roots), intent(out):: weights
        weights = (/1.160440860204388913d-01,&
                                      2.208511247506771413d-01,&
                                      2.413998275878537214d-01,&
                                      1.946367684464170855d-01,&
                                      1.237284159668764899d-01,&
                                      6.367878036898660943d-02,&
                                      2.686047527337972682d-02,&
                                      9.338070881603925677d-03,&
                                      2.680696891336819664d-03,&
                                      6.351291219408556439d-04,&
                                      1.239074599068830081d-04,&
                                      1.982878843895233056d-05,&
                                      2.589350929131392509d-06,&
                                      2.740942840536013206d-07,&
                                      2.332831165025738197d-08,&
                                      1.580745574778327984d-09,&
                                      8.427479123056716393d-11,&
                                      3.485161234907855443d-12,&
                                      1.099018059753451500d-13,&
                                      2.588312664959080167d-15,&
                                      4.437838059840028968d-17,&
                                      5.365918308212045344d-19,&
                                      4.393946892291604451d-21,&
                                      2.311409794388543236d-23,&
                                      7.274588498292248063d-26,&
                                      1.239149701448267877d-28,&
                                      9.832375083105887477d-32,&
                                      2.842323553402700938d-35,&
                                      1.878608031749515392d-39,&
                                      8.745980440465011553d-45/)


    end subroutine gauss_laguerre_30_weights


end module neoclassics_module

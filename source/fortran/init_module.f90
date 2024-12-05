#ifndef INSTALLDIR
#error INSTALLDIR not defined!
#endif

module init_module

   implicit none

contains

   subroutine init_all_module_vars
      !! Initialise all module variables
      !! This is vital to ensure a 'clean' state of Process before a new run starts,
      !! otherwise components of the previous run's state can persist into the new
      !! run. This matters ever since Process is used as a shared library, rather
      !! than a 'run-once' executable.
      use numerics, only: init_numerics
      use process_input, only: init_input
      use buildings_variables, only: init_buildings_variables
      use cost_variables, only: init_cost_variables
      use divertor_variables, only: init_divertor_variables
      use error_handling, only: init_error_handling
      use fson_library, only: init_fson_library
      use fwbs_variables, only: init_fwbs_variables
      use global_variables, only: init_global_variables
      use ccfe_hcpb_module, only: init_ccfe_hcpb_module
      use heat_transport_variables, only: init_heat_transport_variables
      use ife_variables, only: init_ife_variables
      use impurity_radiation_module, only: init_impurity_radiation_module
      use pfcoil_module, only: init_pfcoil_module
      use physics_module, only: init_physics_module
      use physics_variables, only: init_physics_variables
      use scan_module, only: init_scan_module
      use sctfcoil_module, only: init_sctfcoil_module
      use stellarator_module, only: init_stellarator_module
      use stellarator_variables, only: init_stellarator_variables
      use tfcoil_variables, only: init_tfcoil_variables
      use times_variables, only: init_times_variables
      use constants, only: init_constants
      use current_drive_variables, only: init_current_drive_variables
      use primary_pumping_variables, only: init_primary_pumping_variables
      use pfcoil_variables, only: init_pfcoil_variables
      use structure_variables, only: init_structure_variables
      use vacuum_variables, only: init_vacuum_variables
      use pf_power_variables, only: init_pf_power_variables
      use build_variables, only: init_build_variables
      use constraint_variables, only: init_constraint_variables
      use pulse_variables, only: init_pulse_variables
      use rebco_variables, only: init_rebco_variables
      use reinke_variables, only: init_reinke_variables
      use define_iteration_variables, only: init_define_iteration_variables
      use reinke_module, only: init_reinke_module
      use water_usage_variables, only: init_watuse_variables
      use CS_fatigue_variables, only: init_CS_fatigue_variables
      use blanket_library, only: init_blanket_library
      use dcll_module, only: init_dcll_module

      call init_numerics
      call init_input
      call init_buildings_variables
      call init_cost_variables
      call init_divertor_variables
      call init_error_handling
      call init_fson_library
      call init_fwbs_variables
      call init_global_variables
      call init_ccfe_hcpb_module
      call init_heat_transport_variables
      call init_ife_variables
      call init_impurity_radiation_module
      call init_pfcoil_module
      call init_physics_module
      call init_physics_variables
      call init_scan_module
      call init_sctfcoil_module
      call init_stellarator_module
      call init_stellarator_variables
      call init_tfcoil_variables
      call init_times_variables
      call init_constants
      call init_current_drive_variables
      call init_primary_pumping_variables
      call init_pfcoil_variables
      call init_structure_variables
      call init_vacuum_variables
      call init_pf_power_variables
      call init_build_variables
      call init_constraint_variables
      call init_pulse_variables
      call init_rebco_variables
      call init_reinke_variables
      call init_define_iteration_variables
      call init_reinke_module
      call init_watuse_variables
      call init_CS_fatigue_variables
      call init_blanket_library
      call init_dcll_module
   end subroutine init_all_module_vars

   subroutine init

      !! Routine that calls the initialisation routines
      !! author: P J Knight, CCFE, Culham Science Centre
      !! None
      !! This routine calls the main initialisation routines that set
      !! the default values for the global variables, reads in data from
      !! the input file, and checks the run parameters for consistency.

      use global_variables, only: verbose, fileprefix, output_prefix
      use main_module, only: run_summary
      use constants, only: opt_file, vfile, nout, nplot, mfile, sig_file
      use error_handling, only: initialise_error_list
      use numerics, only: ixc , lablxc, nvar
      use process_input, only: nin, input
      use stellarator_module, only: stinit
      implicit none

      !  Arguments

      !  Local variables
      integer :: i

      ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      !  Initialise error handling

      call initialise_error_list

      !  Initialise the program variables
      call initial

      !  Open the input/output external files
      if (trim(fileprefix) == "") then
         open(unit=nin,file="IN.DAT",status='old')
      else
         open(unit=nin,file=trim(fileprefix),status='old')
      end if
      ! open(unit=nin,file=trim(fileprefix)//'IN.DAT',status='old')

      open(unit=nout     ,file=trim(output_prefix)//'OUT.DAT'   ,status='unknown')
      open(unit=mfile    ,file=trim(output_prefix)//'MFILE.DAT' ,status='unknown')

      !  Input any desired new initial values
      call input

      !  Initialise stellarator parameters if necessary
      !  This overrides some of the bounds of the tokamak parameters
      call stinit

      !  Check input data for errors/ambiguities
      call check

      !  Write to the output file certain relevant details about this run
      call run_summary

      !  Open verbose diagnostics file
      if (verbose == 1) then
         open(unit=vfile,file=trim(output_prefix)//'VFILE.DAT',status='unknown')
         write(vfile,'(a80)') 'nviter = number of VMCON iterations.'
         write(vfile,'(a80)') '(1-mod(ifail,7))=1 indicates that there has '// &
            'been an escape from a failed line search.'
         write(vfile,'(a80)') 'odd/even is a convenient plotting bit.'
         write(vfile,'(100a13)') 'nviter','escape', 'odd/even', 'te','coe','rmajor', &
            'fusion_power','bt','t_burn','sqsumsq', (lablxc(ixc(i)),i=1,nvar)
      end if

   end subroutine init

   subroutine open_idempotence_files
      ! Open new output file and mfile to write output to
      ! This is used when checking model evaluation idempotence, to avoid
      ! polluting the final output file and mfile with intermediate result checks
      use global_variables, only: output_prefix
      use constants, only: nout, mfile
      implicit none

      ! Close existing output file and mfile (could be original out and mfiles
      ! or idem scratch files)
      close(unit = nout)
      close(unit = mfile)
      ! Open scratch files with same units
      open(unit=nout, file=trim(output_prefix)//'IDEM_OUT.DAT', action='write', status='replace')
      open(unit=mfile, file=trim(output_prefix)//'IDEM_MFILE.DAT', action='write', status='replace')
   end subroutine open_idempotence_files

   subroutine close_idempotence_files
      ! Close the intermediate idempotence-check files, deleting them in the process
      ! Re-open the original OUT.DAT and MFILE.DAT output files, ready to write
      ! the final data, now model evaluation idempotence has been checked
      use global_variables, only: output_prefix
      use constants, only: nout, mfile
      implicit none

      ! Close idempotence files, deleting them in the process
      close(unit = nout, status="delete")
      close(unit = mfile, status="delete")
      ! Re-open original output file and mfile, appending future output to them
      open(unit=nout, file=trim(output_prefix)//'OUT.DAT', action='write', position='append')
      open(unit=mfile, file=trim(output_prefix)//'MFILE.DAT', action='write', position='append')
   end subroutine close_idempotence_files

   subroutine finish
      ! Originally at the end of the "program", this subroutine writes some final
      ! lines via the output module and then closes any open files. This is
      ! currently called from Python, and will be removed once file handling is
      ! completely dealt with in Python
      ! # TODO Move this output and file handling to Python

      use process_input, only: nin
      use constants, only: iotty, mfile, nout, nplot, opt_file, vfile
      use process_output, only: oheadr
      use global_variables, only: verbose
      implicit none

      call oheadr(nout,'End of PROCESS Output')
      call oheadr(iotty,'End of PROCESS Output')
      call oheadr(nout,'Copy of PROCESS Input Follows')

      close(unit = nin)
      close(unit = nout)
      close(unit = nplot)
      close(unit = mfile)
      close(unit = opt_file)
      if (verbose == 1) close(unit = vfile)
   end subroutine finish

end module init_module

module init_module

#ifndef dp
use, intrinsic :: iso_fortran_env, only: dp=>real64
#endif

   implicit none

   integer, parameter :: nin = 10

contains

   subroutine open_files
      use global_variables, only: verbose, fileprefix, output_prefix
      use constants, only: nout, mfile

      implicit none

      !  Open the input/output external files
      if (trim(fileprefix) == "") then
         open(unit=nin,file="IN.DAT",status='old')
      else
         open(unit=nin,file=trim(fileprefix),status='old')
      end if

      open(unit=nout     ,file=trim(output_prefix)//'OUT.DAT'   ,status='unknown')
      open(unit=mfile    ,file=trim(output_prefix)//'MFILE.DAT' ,status='unknown')

   end subroutine open_files

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

      use constants, only: iotty, mfile, nout, nplot, opt_file, vfile
      use global_variables, only: verbose
      implicit none

      close(unit = nin)
      close(unit = nout)
      close(unit = nplot)
      close(unit = mfile)
      close(unit = opt_file)
      if (verbose == 1) close(unit = vfile)
   end subroutine finish
end module init_module

# Python-Fortran f2py Interface

!!! note
    This documents the abstract Python-Fortran interface which is implemented automatically by f2py upon build. 
    An abstract class, ``foo``, within the abstract interface will be denoted as ``class process._fortran.foo`` while the implementation is accessed using ``process.fortran.foo``; notice the lack of an underscore on the ``fortran`` module.

::: process._fortran
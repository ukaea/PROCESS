# f2py Programming Style Requirements

f2py is a tool developed by Numpy, it is used to provide an interface between Python and Fortran. 
PROCESS uses f2py for its Fortran interface. 

f2py can wrap most modern Fortran, that being said it is opinionated on the structure of modules, 
subroutines and function. 

The following rules will be enforced on all Fortran code, whether currently wrapped or not, and 
any code not adhering to such rules will not pass code review. 

Examples of rules may be given in the context of an example to aide understanding.

## Derived types

f2py does not allow the use of derived types as subroutine/function arguments. For this reason, 
derived types should be avoided or necessary, and instead standard-type arguments passed as parameters.

**Bad:**
```fortran
module foo
    implicit none
    
    type :: duck 
        integer :: id
        character(len=10) :: first_name
    end type duck

    subroutine say_hello(duck_instance)
        type(duck), intent(in) :: duck_instance

        write(*,*) 'Duck ', duck_instance%first_name, ' (with ID ', duck_instance%id, ') says quack!'
    end subroutine say_hello
end module foo
```

**Good:**
```fortran
module foo
    implicit none

    subroutine say_hello(duck_first_name, duck_id)
        character(len=10), intent(in) :: duck_first_name
        integer, intent(in) :: duck_id

        write(*,*) 'Duck ', duck_instance%first_name, ' (with ID ', duck_instance%id, ') says hello!'
    end subroutine say_hello
end module foo
```

## Parameters as dimensions of an array

f2py does not allow dimensions of an array to be parameters without issuing an f2py directive. 

**Bad:**

```fortran
module example
    ! module variables
    integer, parameter :: number = 50

    subroutine example_routine()
        integer, dimension(number, number) :: matrix

        ...
    end subroutine example_routine

end module example
```

**Good:**

```fortran
module example
    ! module variables
    integer, parameter :: number = 50

    subroutine example_routine()
        !f2py integer, intent(aux) :: number
        integer, dimension(number, number) :: matrix

        ...
    end subroutine example_routine

end module example
```

with `!f2py integer, intent(aux) :: number` being the f2py directive.

## Function type decleration

Fortran function method signatures can be defined as follows:

**Bad:**

```fortran
real function example(p1, p2) result(p3)
    ! arguments
    real, intent(in) :: p1, p2
```

which states that `p1` and `p2` are real parameters of the `example` function and `p3` is the real 
result of `example`. This declaration is valid in Fortran, but will be rejected by f2py; f2py cannot 
interpret the type of `p3`. Instead the following patterns should be used:

**Good:**
```fortran
function example(p1, p2) result(p3)
    ! arguments
    real, intent(in) :: p1, p2

    ! return
    real :: p3
```

### Interfaces

Interfaces declaring functions are used throughout PROCESS to allow a function to be an argument of 
another routine (e.g. for integration). f2py requires that function declaration within an interface 
follows the above rules:

**Bad:**

```fortran
interface
    function fun(rho)
        real(dp), intent(in) :: rho
        real(dp) :: fint
    end function fun
end interface
```

Where `fint` is implicitly `intent(out)`. f2py requires an explicit `results` is provided:

**Good:**

```fortran
interface
    function fun(rho) result(fint)
        real(dp), intent(in) :: rho
        real(dp) :: fint
    end function fun
end interface
```

Also note that when using this interface as an argument, the standard Fortran way would be 
`real(dp), external :: fun`. However, when a `results` is given, this is the same as providing 
the `real(dp)` decleration (see above); only `external :: fun` is required.

## Private/public attributes

f2py cannot see `private` routines or variables. This means the use of `private` must be carefully 
used otherwise f2py will raise errors thinking it cannot see variables. This does, however, mean 
that the `private` keyword can be used to 'hide' error-causing variables from f2py - such as 
derived types that are not used as routine parameters. 

Developers should also be carefully not to accidentally implicitly-private a module by using 
the `public` keyword on a variable or routine; every other variable/module variable will become 
implicitly private.

Generally, we advise against the use of `private` and `public` keywords unless in very select circumstances.

## Intrinsic double precision

The Fortran-recommended way to handle `real` precision is to `use` `real64` from the intrinsic `iso_fortran_env`:

**Bad:**

```fortran
! source/fortran/foo.f90
module foo
    use, intrinsic :: iso_fortran_env, only: dp=>real64

    implicit none

    contains 

    subroutine bar(x, y)
        real(dp), intent(in) :: x
        real(dp), intent(out) :: y

        y = x*2
    end subroutine bar
end module foo
```

When f2py wraps this module, it will wrap each subroutine seperately and the `dp` pointer will not 
be available within this wrapper. Because of this, an error will be raised declaring `dp` to be undefined. 

The solution is to preprocess wrapped files first, and replace `dp` with `8`. We use an `ifndef` 
preprocessor directive to indicate not to include the `use ...` line in the wrapping process. 

**Good:**

```fortran
! source/fortran/foo.f90
module foo
#ifndef dp
    use, intrinsic :: iso_fortran_env, only: dp=>real64
#endif

    implicit none

    contains 

    subroutine bar(x, y)
        real(dp), intent(in) :: x
        real(dp), intent(out) :: y

        y = x*2
    end subroutine bar
end module foo
```

This then creates an intermediate source file which does not have undefined pointers:

```fortran
! build/foo.f90
! intermediate source file
module foo
    implicit none

    contains 

    subroutine bar(x, y)
        real(8), intent(in) :: x
        real(8), intent(out) :: y

        y = x*2
    end subroutine bar
end module foo
```

## Assumed-shape arrays

f2py needs to have a definite size for each dimension of an array to be passed back through the 
Python-Fortran interface, ie `intent(out)`. This means assumed-shape and assumed-size arrays will 
not work when declared as a return value.

**Bad:**

```fortran
module foo
    use, intrinsic :: iso_fortran_env, only: dp=>real64

    implicit none

    contains 

    subroutine bar(x, y)
        real(dp), dimension(:,:), intent(in) :: x
        real(dp), dimension(:,:), intent(out) :: y

        ! routine logic
    end subroutine bar
end module foo
```

For example, if both arrays have dimensions of sizes `n` and `i`, respectively.

**Good:**

```fortran
module foo
    use, intrinsic :: iso_fortran_env, only: dp=>real64

    implicit none

    contains 

    subroutine bar(x, y, n, i)
        integer, intent(in) :: n, i
        real(dp), dimension(:,:), intent(in) :: x
        real(dp), dimension(n,i), intent(out) :: y

        ! routine logic
    end subroutine bar
end module foo
```

Notice that the input array, `x`, does not need to be changed.

## Strings

### Length

Firstly, a string being set via the Python-Fortran interface must be of the exact length specified, 
f2py will reject anything under or over. For this reason, we provide a utility function 
`process.utilities.f2py_string_patch.string_to_f2py_compatible`, to be used as follows:

```python
! foo.py

from process.utilities.f2py_string_patch import string_to_f2py_compatible


def bar(string: str) -> None:
    fortran.foomod.barvar = string_to_f2py_compatible(fortran.foomod.barvar, string)
```

where `fortran.foomod.barvar` will be set as the string `string`, padded with whitespace if too 
short or truncated if too long. 

Remember to also use the fortran `trim` function to remove padding when using any strings inside 
of Fortran code.

### Arrays

Character arrays in Fortran, declared like `character(len=10), dimension(m), intent(out) :: foo`, 
are not liked by f2py. Instead, the slightly outdated `character*10, dimension(m), intent(out) :: foo` 
syntax is needed. This error will manifest itself as the character array (`foo`) only being of size 
1 regardless of the value of `m`.

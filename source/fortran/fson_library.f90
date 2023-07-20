! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!+PJK Need to convert to autodoc style + usual PROCESS standard layout

!-----------------------------------------------------------------------------------------
!
!  FSON library
!
!  Extracted from https://github.com/josephalevin/fson on 9th July 2014
!
!    with selected updates taken from the version forked as
!
!  https://github.com/jmozmoz/fson/commit/b210a9011bb804957546e2a4b6eade578e7035ef
!
!    plus some improvements to help with array handling and double precision, by
!    P J Knight, 17th July 2014
!
! Comprises the following original FSON files:
!    string.f95, value_m.f95, fson_path_m.f95, fson.f95
!
!-----------------------------------------------------------------------------------------
!
! Copyright (c) 2012 Joseph A. Levin
!
! Permission is hereby granted, free of charge, to any person obtaining a copy of this
! software and associated documentation files (the "Software"), to deal in the Software
! without restriction, including without limitation the rights to use, copy, modify, merge,
! publish, distribute, sublicense, and/or sell copies of the Software, and to permit
! persons to whom the Software is furnished to do so, subject to the following conditions:
!
! The above copyright notice and this permission notice shall be included in all copies or
! substantial portions of the Software.
!
! THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED,
! INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR
! PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE
! LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT
! OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
! DEALINGS IN THE SOFTWARE.
!
!-----------------------------------------------------------------------------------------

module fson_string_m

  private

  public :: fson_string, fson_string_create, fson_string_destroy, fson_string_length, &
       fson_string_append, fson_string_clear
  public :: fson_string_equals, fson_string_copy

  integer, parameter :: BLOCK_SIZE = 32

  type fson_string
     character (len = BLOCK_SIZE) :: chars
     integer :: index = 0
     type(fson_string), pointer :: next => null()
  end type fson_string

  interface fson_string_append
     module procedure append_chars, append_string
  end interface

  interface fson_string_copy
     module procedure copy_chars
  end interface

  interface fson_string_equals
     module procedure equals_string
  end interface

  interface fson_string_length
     module procedure string_length
  end interface

contains

  !
  ! FSON STRING CREATE
  !
  function fson_string_create(chars) result(new)
    character(len=*), optional :: chars
    type(fson_string), pointer :: new

    allocate(new)

    ! append chars if available
    if (present(chars)) then
       call append_chars(new, chars)
    end if

  end function fson_string_create

  !
  ! FSON STRING DESTROY
  !
  recursive subroutine fson_string_destroy(this)
    type(fson_string), pointer :: this

    if (associated(this % next)) then
       call fson_string_destroy(this % next)
    end if

    nullify (this % next)
    nullify (this)
  end subroutine fson_string_destroy

  !
  ! ALLOCATE BLOCK
  !
  subroutine allocate_block(this)
    type(fson_string), pointer :: this
    type(fson_string), pointer :: new

    if (.not.associated(this % next)) then
       allocate(new)
       this % next => new
    end if

  end subroutine allocate_block

  !
  ! APPEND_STRING
  !
  subroutine append_string(str1, str2)
    type(fson_string), pointer :: str1, str2
    integer length, i

    length = string_length(str2)

    do i = 1, length
       call append_char(str1, get_char_at(str2, i))
    end do

  end subroutine append_string

  !
  ! APPEND_CHARS
  !
  subroutine append_chars(str, c)
    type(fson_string), pointer :: str
    character (len = *), intent(in) :: c
    integer length, i

    length = len(c)

    do i = 1, length
       call append_char(str, c(i:i))
    end do

  end subroutine append_chars

  !
  ! APPEND_CHAR
  !
  recursive subroutine append_char(str, c)
    type(fson_string), pointer :: str
    character, intent(in) :: c

    if (str % index >= BLOCK_SIZE) then
       !set down the chain
       call allocate_block(str)
       call append_char(str % next, c)

    else
       ! set local
       str % index = str % index + 1
       str % chars(str % index:str % index) = c
    end if

  end subroutine append_char

  !
  ! COPY CHARS
  !
  subroutine copy_chars(this, to)
    type(fson_string), pointer :: this
    character(len = *), intent(inout) :: to
    integer :: length

    length = min(string_length(this), len(to))

    do i = 1, length
       to(i:i) = get_char_at(this, i)
    end do

    ! pad with nothing
    do i = length + 1, len(to)
       to(i:i) = ""
    end do

  end subroutine copy_chars

  !
  ! CLEAR
  !
  recursive subroutine string_clear(this)
    type(fson_string), pointer :: this

    if (associated(this % next)) then
       call string_clear(this % next)
       deallocate(this % next)
       nullify (this % next)
    end if

    this % index = 0

  end subroutine string_clear

  !
  ! SIZE
  !
  recursive integer function string_length(str) result(count)
    type(fson_string), pointer :: str

    count = str % index

    if (str % index == BLOCK_SIZE .AND. associated(str % next)) then
       count = count + string_length(str % next)
    end if

  end function string_length

  !
  ! GET CHAR AT
  !
  recursive character function get_char_at(this, i) result(c)
    type(fson_string), pointer :: this
    integer, intent(in) :: i

    if (i <= this % index) then
       c = this % chars(i:i)
    else
       c = get_char_at(this % next, i - this % index)
    end if

  end function get_char_at

  !
  ! EQUALS STRING
  !
  logical function equals_string(this, other) result(equals)
    type(fson_string), pointer :: this, other
    integer :: i
    equals = .false.

    if (fson_string_length(this) /= fson_string_length(other)) then
       equals = .false.
       return
    else if (fson_string_length(this) == 0) then
       equals = .true.
       return
    end if

    do i=1, fson_string_length(this)
       if (get_char_at(this, i) /= get_char_at(other, i)) then
          equals = .false.
          return
       end if
    end do

    equals = .true.

  end function equals_string

end module fson_string_m

! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module fson_value_m
  use fson_string_m, only: fson_string
  implicit none

  private

  public :: fson_value, fson_value_create, fson_value_destroy, fson_value_add, fson_value_get, &
       fson_value_count, fson_value_print

  ! constants for the value types
  integer, public, parameter :: TYPE_UNKNOWN = -1
  integer, public, parameter :: TYPE_NULL = 0
  integer, public, parameter :: TYPE_OBJECT = 1
  integer, public, parameter :: TYPE_ARRAY = 2
  integer, public, parameter :: TYPE_STRING = 3
  integer, public, parameter :: TYPE_INTEGER = 4
  integer, public, parameter :: TYPE_REAL = 5
  integer, public, parameter :: TYPE_LOGICAL = 6

  !
  ! FSON VALUE
  !
  type fson_value
     type(fson_string), pointer :: name => null()
     integer :: value_type = TYPE_UNKNOWN
     logical :: value_logical
     integer :: value_integer
     real :: value_real
     !+PJK
     real(kind(1.0D0)) :: value_double
     !-PJK
     type(fson_string), pointer :: value_string => null()
     type(fson_value), pointer :: next => null()
     type(fson_value), pointer :: parent => null()
     type(fson_value), pointer :: children => null()
  end type fson_value

  !
  ! FSON VALUE GET
  !
  ! Use either a 1 based index or member name to get the value.
  interface fson_value_get
     module procedure get_by_index
     module procedure get_by_name_chars
     module procedure get_by_name_string
  end interface

contains

  !
  ! FSON VALUE CREATE
  !
  function fson_value_create() result(new)
    type(fson_value), pointer :: new

    allocate(new)

  end function fson_value_create

  !
  ! FSON VALUE DESTROY
  !
  recursive subroutine fson_value_destroy(this)
    use fson_string_m, only: fson_string_destroy
    implicit none

    type(fson_value), pointer :: this

    if (associated(this % children)) then
       call fson_value_destroy(this % children)
       nullify(this % children)
    end if

    if (associated(this % next)) then
       call fson_value_destroy(this % next)
       nullify (this % next)
    end if

    if (associated(this % name)) then
       call fson_string_destroy(this % name)
       nullify (this % name)
    end if

    if (associated(this % value_string)) then
       call fson_string_destroy(this % value_string)
       nullify (this % value_string)
    end if

    nullify(this)

  end subroutine fson_value_destroy

  !
  ! FSON VALUE ADD
  !
  ! Adds the member to the linked list
  subroutine fson_value_add(this, member)
    type(fson_value), pointer :: this, member, p

    ! associate the parent
    member % parent => this

    ! add to linked list
    if (associated(this % children)) then
       ! get to the tail of the linked list
       p => this % children
       do while (associated(p % next))
          p => p % next
       end do

       p % next => member
    else
       this % children => member
    end if

  end subroutine fson_value_add

  !
  ! FSON_VALUE_COUNT
  !
  integer function fson_value_count(this) result(count)
    type(fson_value), pointer :: this, p

    count = 0

    p => this % children

    do while (associated(p))
       count = count + 1
       p => p % next
    end do

  end function fson_value_count

  !
  ! GET BY INDEX
  !
  function get_by_index(this, index) result(p)
    type(fson_value), pointer :: this, p
    integer, intent(in) :: index
    integer :: i

    p => this % children

    do i = 1, index - 1
       p => p % next
    end do

  end function get_by_index

  !
  ! GET BY NAME CHARS
  !
  function get_by_name_chars(this, name) result(p)
    use fson_string_m, only: fson_string, fson_string_create
    implicit none

    type(fson_value), pointer :: this, p
    character(len=*), intent(in) :: name

    type(fson_string), pointer :: string

    ! convert the char array into a string
    string => fson_string_create(name)

    p => get_by_name_string(this, string)

  end function get_by_name_chars

  !
  ! GET BY NAME STRING
  !
  function get_by_name_string(this, name) result(p)
    use fson_string_m, only: fson_string, fson_string_equals
    implicit none

    type(fson_value), pointer :: this, p
    type(fson_string), pointer :: name
    integer :: i

    if (this % value_type /= TYPE_OBJECT) then
       nullify(p)
       return
    end if

    do i=1, fson_value_count(this)
       p => fson_value_get(this, i)
       if (fson_string_equals(p%name, name)) then
          return
       end if
    end do

    ! didn't find anything
    nullify(p)

  end function get_by_name_string

  !
  ! FSON VALUE PRINT
  !
  recursive subroutine fson_value_print(this, indent)
    use fson_string_m, only: fson_string_copy
    implicit none

    type(fson_value), pointer :: this, element
    integer, optional, intent(in) :: indent
    character(len=1024) :: tmp_chars
    integer :: tab, i, count, spaces

    !+PJK
    if (.not.associated(this)) return
    !-PJK

    if (present(indent)) then
       tab = indent
    else
       tab = 0
    end if

    spaces = tab * 2

    select case (this % value_type)
    case(TYPE_OBJECT)
       print *, repeat(" ", spaces), "{"
       count = fson_value_count(this)
       do i = 1, count
          ! get the element
          element => fson_value_get(this, i)
          ! get the name
          call fson_string_copy(element % name, tmp_chars)
          ! print the name
          print *, repeat(" ", spaces), '"', trim(tmp_chars), '":'
          ! recursive print of the element
          call fson_value_print(element, tab + 1)
          ! print the separator if required
          if (i < count) then
             print *, repeat(" ", spaces), ","
          end if
       end do

       print *, repeat(" ", spaces), "}"
    case (TYPE_ARRAY)
       print *, repeat(" ", spaces), "["
       count = fson_value_count(this)
       do i = 1, count
          ! get the element
          element => fson_value_get(this, i)
          ! recursive print of the element
          call fson_value_print(element, tab + 1)
          ! print the separator if required
          if (i < count) then
             print *, ","
          end if
       end do
       print *, repeat(" ", spaces), "]"
    case (TYPE_NULL)
       print *, repeat(" ", spaces), "null"
    case (TYPE_STRING)
       call fson_string_copy(this % value_string, tmp_chars)
       print *, repeat(" ", spaces), '"', trim(tmp_chars), '"'
    case (TYPE_LOGICAL)
       if (this % value_logical) then
          print *, repeat(" ", spaces), "true"
       else
          print *, repeat(" ", spaces), "false"
       end if
    case (TYPE_INTEGER)
       print *, repeat(" ", spaces), this % value_integer
    case (TYPE_REAL)
       print *, repeat(" ", spaces), this % value_real  !  N.B. doubles will be shown as single precision
    end select
  end subroutine fson_value_print

end module fson_value_m

! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Modifications by P J Knight to handle arrays more conveniently than by providing
! an array_callback argument to get_array, and to handle double precision values better

module fson_path_m

  private

  public :: fson_path_get, array_callback

  interface fson_path_get
     module procedure get_by_path
     module procedure get_integer
     module procedure get_real
     module procedure get_double
     module procedure get_logical
     module procedure get_chars
     module procedure get_array
     !+PJK
     module procedure get_int_array
     module procedure get_real_array
     module procedure get_double_array
     module procedure get_string_array
     module procedure get_int_array_in_struct
     module procedure get_real_array_in_struct
     module procedure get_double_array_in_struct
     module procedure get_string_array_in_struct
     !-PJK
  end interface

contains

  !
  ! GET BY PATH
  !
  ! $     = root
  ! @     = this
  ! .     = child object member
  ! []    = child array element
  !
  recursive subroutine get_by_path(this, path, p)
    use fson_value_m, only: fson_value, fson_value_get
    implicit none

    type(fson_value), pointer :: this, p
    character(len=*), intent(inout) :: path
    integer :: i, length, child_i
    character :: c
    logical :: array

    ! default to assuming relative to this
    p => this

    child_i = 1

    array = .false.

    length = len_trim(path)

    do i=1, length
       c = path(i:i)
       select case (c)
       case ("$")
          ! root
          do while (associated (p % parent))
             p => p % parent
          end do
          child_i = i + 1
       case ("@")
          ! this
          p => this
          child_i = i + 1
       case (".", "[")
          ! get child member from p
          if (child_i < i) then
             p => fson_value_get(p, path(child_i:i-1))
          else
             child_i = i + 1
             cycle
          end if

          if (.not.associated(p)) then
             return
          end if

          child_i = i + 1

          ! check if this is an array
          ! if so set the array flag
          if (c == "[") then
             ! start looking for the array element index
             array = .true.
          end if
       case ("]")
          if (.not.array) then
             print *, "ERROR: Unexpected ], not missing preceding ["
             call exit(1)
          end if
          array = .false.
          child_i = parse_integer(path(child_i:i-1))
          p => fson_value_get(p, child_i)

          child_i = i + 1
       end select
    end do

    ! grab the last child if present in the path
    if (child_i <= length) then
       p => fson_value_get(p, path(child_i:i-1))
       if (.not.associated(p)) then
          return
       else
       end if
    end if

  end subroutine get_by_path

  !
  ! PARSE INTEGER
  !
  integer function parse_integer(chars) result(integral)
    character(len=*) :: chars
    character :: c
    integer :: tmp, i

    integral = 0
    do i=1, len_trim(chars)
       c = chars(i:i)
       select case(c)
       case ("0":"9")
          ! digit
          read (c, '(i1)') tmp
          ! shift
          if (i > 1) then
             integral = integral * 10
          end if
          ! add
          integral = integral + tmp

       case default
          return
       end select
    end do

  end function parse_integer

  !
  ! GET INTEGER
  !
  subroutine get_integer(this, path, value)
    use fson_value_m, only: type_integer, type_real, type_logical, fson_value
    implicit none

    type(fson_value), pointer :: this, p
    character(len=*), optional :: path
    integer :: value

    nullify(p)
    if (present(path)) then
       call get_by_path(this=this, path=path, p=p)
    else
       p => this
    end if

    if (.not.associated(p)) then
       print *, "Unable to resolve path: ", path
       call exit(1)
    end if

    if (p % value_type == TYPE_INTEGER) then
       value = p % value_integer
    else if (p % value_type == TYPE_REAL) then
       value = p % value_real
    else if (p % value_type == TYPE_LOGICAL) then
       if (p % value_logical) then
          value = 1
       else
          value = 0
       end if
    else
       print *, "Unable to resolve value to integer: ", path
       call exit(1)
    end if

  end subroutine get_integer

  !
  ! GET REAL
  !
  subroutine get_real(this, path, value)
    use fson_value_m, only: type_integer, type_real, type_logical, fson_value
    implicit none

    type(fson_value), pointer :: this, p
    character(len=*), optional :: path
    real :: value

    nullify(p)

    if (present(path)) then
       call get_by_path(this=this, path=path, p=p)
    else
       p => this
    end if

    if (.not.associated(p)) then
       print *, "Unable to resolve path: ", path
       call exit(1)
    end if

    if (p % value_type == TYPE_INTEGER) then
       value = p % value_integer
    else if (p % value_type == TYPE_REAL) then
       value = p % value_real
    else if (p % value_type == TYPE_LOGICAL) then
       if (p % value_logical) then
          value = 1
       else
          value = 0
       end if
    else
       print *, "Unable to resolve value to real: ", path
       call exit(1)
    end if

  end subroutine get_real

  !
  ! GET DOUBLE
  !
  subroutine get_double(this, path, value)
    use fson_value_m, only: type_integer, type_real, type_logical, fson_value
    implicit none

    type(fson_value), pointer :: this, p
    character(len=*), optional :: path
    real(kind(1.0D0)) :: value

    nullify(p)

    if (present(path)) then
       call get_by_path(this=this, path=path, p=p)
    else
       p => this
    end if

    if (.not.associated(p)) then
       print *, "Unable to resolve path: ", path
       call exit(1)
    end if

    if (p % value_type == TYPE_INTEGER) then
       value = p % value_integer
    else if (p % value_type == TYPE_REAL) then
       value = p % value_double  !  PJK from value_real
    else if (p % value_type == TYPE_LOGICAL) then
       if (p % value_logical) then
          value = 1
       else
          value = 0
       end if
    else
       print *, "Unable to resolve value to double: ", path
       call exit(1)
    end if

  end subroutine get_double

  !
  ! GET LOGICAL
  !
  subroutine get_logical(this, path, value)
    use fson_value_m, only: type_integer, type_logical, fson_value
    implicit none

    type(fson_value), pointer :: this, p
    character(len=*), optional :: path
    logical :: value

    nullify(p)

    if (present(path)) then
       call get_by_path(this=this, path=path, p=p)
    else
       p => this
    end if

    if (.not.associated(p)) then
       print *, "Unable to resolve path: ", path
       call exit(1)
    end if

    if (p % value_type == TYPE_INTEGER) then
       value = (p % value_integer > 0)
    else if (p % value_type == TYPE_LOGICAL) then
       value = p % value_logical
    else
       print *, "Unable to resolve value to logical: ", path
       call exit(1)
    end if

  end subroutine get_logical

  !
  ! GET CHARS
  !
  subroutine get_chars(this, path, value)
    use fson_value_m, only: type_string, fson_value
    use fson_string_m, only: fson_string_copy
    implicit none

    type(fson_value), pointer :: this, p
    character(len=*), optional :: path
    character(len=*) :: value

    nullify(p)

    if (present(path)) then
       call get_by_path(this=this, path=path, p=p)
    else
       p => this
    end if

    if (.not.associated(p)) then
       print *, "Unable to resolve path: ", path
       call exit(1)
    end if

    if (p % value_type == TYPE_STRING) then
       call fson_string_copy(p % value_string, value)
    else
       print *, "Unable to resolve value to characters: ", path
       call exit(1)
    end if

  end subroutine get_chars

  !
  ! GET ARRAY (original version using array_callback)
  !
  subroutine get_array(this, path, array_callback)
    use fson_value_m, only: type_array, fson_value_get, fson_value_count, &
      fson_value
    implicit none

    type(fson_value), pointer :: this, p, element
    character(len=*), optional :: path
    integer :: index, count

    ! ELEMENT CALLBACK  (PJK: Added example comments)
    interface
       subroutine array_callback(element, index, count)
         use fson_value_m, only: fson_value
         implicit none

         !  In the actual routine add a second 'use' line as follows:
         !use shared_data  !  contains declarations for the array(s) to be populated

         type(fson_value), pointer :: element
         integer :: index, count

         !  Example usage in the actual routine:
         !  call fson_get(element, "element_name", array_name(index))

       end subroutine array_callback
    end interface

    nullify(p)

    ! resolve the path to the value
    if (present(path)) then
       call get_by_path(this=this, path=path, p=p)
    else
       p => this
    end if

    if (.not.associated(p)) then
       print *, "Unable to resolve path: ", path
       call exit(1)
    end if

    if (p % value_type == TYPE_ARRAY) then
       count = fson_value_count(p)
       do index = 1, count
          element => fson_value_get(p, index)
          call array_callback(element, index, count)
       end do
    else
       print *, "Resolved value is not an array. ", path
       call exit(1)
    end if

  end subroutine get_array

  !+PJK
  !
  ! GET INT ARRAY
  !
  subroutine get_int_array(this, path, array)
    use fson_value_m, only: type_array, fson_value_get, fson_value_count, &
      fson_value
    implicit none

    type(fson_value), pointer :: this, p, element
    character(len=*), optional :: path
    integer :: index, count
    integer, dimension(:) :: array

    nullify(p)

    ! resolve the path to the value
    if (present(path)) then
       call get_by_path(this=this, path=path, p=p)
    else
       p => this
    end if

    if (.not.associated(p)) then
       print *, "Unable to resolve path: ", path
       call exit(1)
    end if

    if (p % value_type == TYPE_ARRAY) then
       count = fson_value_count(p)
       do index = 1, count
          element => fson_value_get(p, index)
          array(index) = element%value_integer
       end do
    else
       print *, "Resolved value is not an array. ", path
       call exit(1)
    end if

  end subroutine get_int_array

  !
  ! GET REAL ARRAY
  !
  subroutine get_real_array(this, path, array)
    use fson_value_m, only: type_array, fson_value_get, fson_value_count, &
      fson_value
    implicit none

    type(fson_value), pointer :: this, p, element
    character(len=*), optional :: path
    integer :: index, count
    real, dimension(:) :: array

    nullify(p)

    ! resolve the path to the value
    if (present(path)) then
       call get_by_path(this=this, path=path, p=p)
    else
       p => this
    end if

    if (.not.associated(p)) then
       print *, "Unable to resolve path: ", path
       call exit(1)
    end if

    if (p % value_type == TYPE_ARRAY) then
       count = fson_value_count(p)
       do index = 1, count
          element => fson_value_get(p, index)
          array(index) = element%value_real
       end do
    else
       print *, "Resolved value is not an array. ", path
       call exit(1)
    end if

  end subroutine get_real_array

  !
  ! GET DOUBLE ARRAY
  !
  subroutine get_double_array(this, path, array)
    use fson_value_m, only: type_array, fson_value_get, fson_value_count, &
      fson_value
    implicit none

    type(fson_value), pointer :: this, p, element
    character(len=*), optional :: path
    integer :: index, count
    real(kind(1.0D0)), dimension(:) :: array

    nullify(p)

    ! resolve the path to the value
    if (present(path)) then
       call get_by_path(this=this, path=path, p=p)
    else
       p => this
    end if

    if (.not.associated(p)) then
       print *, "Unable to resolve path: ", path
       call exit(1)
    end if

    if (p % value_type == TYPE_ARRAY) then
       count = fson_value_count(p)
       do index = 1, count
          element => fson_value_get(p, index)
          array(index) = element%value_double
       end do
    else
       print *, "Resolved value is not an array. ", path
       call exit(1)
    end if

  end subroutine get_double_array

  !
  ! GET STRING ARRAY
  !
  subroutine get_string_array(this, path, array)
    use fson_value_m, only: fson_value_count, fson_value_get, fson_value, &
      TYPE_ARRAY
    use fson_string_m, only: fson_string_copy
    implicit none

    type(fson_value), pointer :: this, p, element
    character(len=*), optional :: path
    integer :: index, count
    character(len=*), dimension(:) :: array

    nullify(p)

    ! resolve the path to the value
    if (present(path)) then
       call get_by_path(this=this, path=path, p=p)
    else
       p => this
    end if

    if (.not.associated(p)) then
       print *, "Unable to resolve path: ", path
       call exit(1)
    end if

    if (p % value_type == TYPE_ARRAY) then
       count = fson_value_count(p)
       do index = 1, count
          element => fson_value_get(p, index)
          call fson_string_copy(element%value_string, array(index))
       end do
    else
       print *, "Resolved value is not an array. ", path
       call exit(1)
    end if

  end subroutine get_string_array

  !
  ! GET INT ARRAY IN STRUCTURE
  !
  subroutine get_int_array_in_struct(this, path, subpath, array)
    use fson_value_m, only: fson_value_count, fson_value_get, fson_value, &
      TYPE_ARRAY
    implicit none

    type(fson_value), pointer :: this, p, element
    character(len=*) :: path, subpath
    integer, dimension(:), intent(out) :: array
    integer :: index, count

    nullify(p)

    ! resolve the path to the value
    call get_by_path(this=this, path=path, p=p)

    if (.not.associated(p)) then
       print *, "Unable to resolve path: ", path
       call exit(1)
    end if

    if (p % value_type == TYPE_ARRAY) then
       count = fson_value_count(p)
       do index = 1, count
          element => fson_value_get(p, index)
          call get_integer(element, subpath, array(index))
       end do
    else
       print *, "Resolved value is not an array. ", path
       call exit(1)
    end if

  end subroutine get_int_array_in_struct

  !
  ! GET REAL ARRAY IN STRUCTURE
  !
  subroutine get_real_array_in_struct(this, path, subpath, array)
    use fson_value_m, only: fson_value_count, fson_value_get, fson_value, &
      TYPE_ARRAY
    implicit none

    type(fson_value), pointer :: this, p, element
    character(len=*) :: path, subpath
    real, dimension(:), intent(out) :: array
    integer :: index, count

    nullify(p)

    ! resolve the path to the value
    call get_by_path(this=this, path=path, p=p)

    if (.not.associated(p)) then
       print *, "Unable to resolve path: ", path
       call exit(1)
    end if

    if (p % value_type == TYPE_ARRAY) then
       count = fson_value_count(p)
       do index = 1, count
          element => fson_value_get(p, index)
          call get_real(element, subpath, array(index))
       end do
    else
       print *, "Resolved value is not an array. ", path
       call exit(1)
    end if

  end subroutine get_real_array_in_struct

  !
  ! GET DOUBLE ARRAY IN STRUCTURE
  !
  subroutine get_double_array_in_struct(this, path, subpath, array)
    use fson_value_m, only: fson_value_count, fson_value_get, fson_value, &
      TYPE_ARRAY
    implicit none

    type(fson_value), pointer :: this, p, element
    character(len=*) :: path, subpath
    real(kind(1.0D0)), dimension(:), intent(out) :: array
    integer :: index, count

    nullify(p)

    ! resolve the path to the value
    call get_by_path(this=this, path=path, p=p)

    if (.not.associated(p)) then
       print *, "Unable to resolve path: ", path
       call exit(1)
    end if

    if (p % value_type == TYPE_ARRAY) then
       count = fson_value_count(p)
       do index = 1, count
          element => fson_value_get(p, index)
          call get_double(element, subpath, array(index))
       end do
    else
       print *, "Resolved value is not an array. ", path
       call exit(1)
    end if

  end subroutine get_double_array_in_struct

  !
  ! GET STRING ARRAY IN STRUCTURE
  !
  subroutine get_string_array_in_struct(this, path, subpath, array)
    use fson_value_m, only: fson_value_count, fson_value_get, fson_value, &
      TYPE_ARRAY
    implicit none

    type(fson_value), pointer :: this, p, element
    character(len=*) :: path, subpath
    character(len=*), dimension(:), intent(out) :: array
    integer :: index, count

    nullify(p)

    ! resolve the path to the value
    call get_by_path(this=this, path=path, p=p)

    if (.not.associated(p)) then
       print *, "Unable to resolve path: ", path
       call exit(1)
    end if

    if (p % value_type == TYPE_ARRAY) then
       count = fson_value_count(p)
       do index = 1, count
          element => fson_value_get(p, index)
          call get_chars(element, subpath, array(index))
       end do
    else
       print *, "Resolved value is not an array. ", path
       call exit(1)
    end if

  end subroutine get_string_array_in_struct
  !-PJK

end module fson_path_m

! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module fson_library

  !! JSON file reading module
  !! author: P J Knight, CCFE, Culham Science Centre
  !! N/A
  !! This module uses a local copy of the freely-available FSON library,
  !! written by Joseph A. Levin and distributed via github,
  !! to enable PROCESS to read in information from JSON-formatted files.
  !! None
  !
  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  use fson_value_m, only: fson_print => fson_value_print, &
    fson_destroy => fson_value_destroy, fson_value
  use fson_path_m, only: fson_get => fson_path_get

  implicit none

  private

  public :: fson_parse, fson_value, fson_get, fson_print, fson_destroy, &
    init_fson_library

  ! FILE IOSTAT CODES
  integer, parameter :: end_of_file = -1
  integer, parameter :: end_of_record = -2

  ! PARSING STATES
  integer, parameter :: STATE_LOOKING_FOR_VALUE = 1
  integer, parameter :: STATE_IN_OBJECT = 2
  integer, parameter :: STATE_IN_PAIR_NAME = 3
  integer, parameter :: STATE_IN_PAIR_VALUE = 4

  ! POP/PUSH CHARACTER
  integer :: pushed_index
  character (len = 10) :: pushed_char

contains

  subroutine init_fson_library
    !! Initialise fson library module variables
    implicit none

    pushed_index = 0
  end subroutine init_fson_library

  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !
  ! FSON PARSE
  !
  function fson_parse(file, unit) result(p)
    use fson_value_m, only: fson_value_create
    implicit none

    type(fson_value), pointer :: p
    integer, optional, intent(inout) :: unit
    character(len = *), intent(in) :: file
    logical :: unit_available
    integer :: u
    ! init the pointer to null
    nullify(p)

    ! select the file unit to use
    if (present(unit)) then
       u = unit
    else
       ! find the first available unit
       unit_available = .false.
       u = 20

       do while (.not.unit_available)
          inquire(unit=u, exist=unit_available)
          u = u + 1
       end do
    end if

    ! open the file
    open (unit=u, file=file, status="old", action="read", &
         form="formatted", position="rewind")

    ! create the value and associate the pointer
    p => fson_value_create()

    ! parse as a value
    call parse_value(unit=u, value=p)

    ! close the file
    if ( .not. present(unit)) then
       close (u)
    end if

  end function fson_parse

  !
  ! PARSE_VALUE
  !
  recursive subroutine parse_value(unit, value)
    use fson_value_m, only: TYPE_ARRAY, TYPE_LOGICAL, TYPE_NULL, TYPE_OBJECT, &
      TYPE_STRING
    implicit none

    integer, intent(inout) :: unit
    type(fson_value), pointer :: value
    logical :: eof
    character :: c

    ! for some unknown reason the next pointer is getting messed with the pop
    type(fson_value), pointer :: hack

    ! start the hack
    hack => value % next

    ! pop the next non whitespace character off the file
    c = pop_char(unit, eof=eof, skip_ws=.true.)

    ! finish the hack; set the next pointer to whatever it was before the pop
    value % next => hack

    if (eof) then
       return
    else
       select case (c)
       case ("{")
          ! start object
          value % value_type = TYPE_OBJECT
          call parse_object(unit, value)
       case ("[")
          ! start array
          value % value_type = TYPE_ARRAY
          call parse_array(unit, value)
       case ("]")
          ! end an empty array
          nullify(value)
       case ('"')
          ! string
          value % value_type = TYPE_STRING
          value % value_string => parse_string(unit)
       case ("t")
          ! true
          value % value_type = TYPE_LOGICAL
          call parse_for_chars(unit, "rue")
          value % value_logical = .true.
       case ("f")
          ! false
          value % value_type = TYPE_LOGICAL
          value % value_logical = .false.
          call parse_for_chars(unit, "alse")
       case ("n")
          value % value_type = TYPE_NULL
          call parse_for_chars(unit, "ull")
       case("-", "0": "9")
          call push_char(c)
          call parse_number(unit, value)
       case default
          print *, "ERROR: Unexpected character while parsing value. '", c, "' ASCII=", iachar(c)
          call exit (1)
       end select
    end if

  end subroutine parse_value

  !
  ! PARSE OBJECT
  !
  recursive subroutine parse_object(unit, parent)
  use fson_value_m, only: fson_value_create, fson_value_add
  implicit none

    integer, intent(inout) :: unit
    type(fson_value), pointer :: parent, pair

    logical :: eof
    character :: c

    ! pair name
    c = pop_char(unit, eof=eof, skip_ws=.true.)
    if (eof) then
       print *, "ERROR: Unexpected end of file while parsing start of object."
       call exit (1)
    else if ("}" == c) then
       ! end of an empty object
       return
    else if ('"' == c) then
       pair => fson_value_create()
       pair % name => parse_string(unit)
    else
       print *, "ERROR: Expecting string: '", c, "'"
       call exit (1)
    end if

    ! pair value
    c = pop_char(unit, eof=eof, skip_ws=.true.)
    if (eof) then
       print *, "ERROR: Unexpected end of file while parsing object member. 1"
       call exit (1)
    else if (":" == c) then
       ! parse the value
       call parse_value(unit, pair)
       call fson_value_add(parent, pair)
    else
       print *, "ERROR: Expecting : and then a value. ", c
       call exit (1)
    end if

    ! another possible pair
    c = pop_char(unit, eof=eof, skip_ws=.true.)
    if (eof) then
       return
    else if ("," == c) then
       ! read the next member
       call parse_object(unit=unit, parent=parent)
    else if ("}" == c) then
       return
    else
       print *, "ERROR: Expecting end of object.", c
       call exit (1)
    end if

  end subroutine parse_object

  !
  ! PARSE ARRAY
  !
  recursive subroutine parse_array(unit, array)
    use fson_value_m, only: fson_value_create, fson_value_add
    implicit none

    integer, intent(inout) :: unit
    type(fson_value), pointer :: array, element

    logical :: eof
    character :: c

    ! try to parse an element value
    element => fson_value_create()
    call parse_value(unit, element)

    ! parse value will disassociate an empty array value
    if (associated(element)) then
       call fson_value_add(array, element)
    end if

    ! popped the next character
    c = pop_char(unit, eof=eof, skip_ws=.true.)

    if (eof) then
       return
    else if ("," == c) then
       ! parse the next element
       call parse_array(unit, array)
    else if ("]" == c) then
       ! end of array
       return
    end if

  end subroutine parse_array

  !
  ! PARSE STRING
  !
  function parse_string(unit) result(string)
    use fson_string_m, only: fson_string, fson_string_create, fson_string_append
    implicit none

    integer, intent(inout) :: unit
    type(fson_string), pointer :: string

    logical :: eof
    character :: c, last

    string => fson_string_create()

    do
       c = pop_char(unit, eof=eof, skip_ws=.false.)
       if (eof) then
          print *, "Expecting end of string"
          call exit(1)!
       else if ('"' == c .and. last /= '\') then  !'
          exit
       else
          last = c
          call fson_string_append(string, c)
       end if
    end do
  end function parse_string

  !
  ! PARSE FOR CHARACTERS
  !
  subroutine parse_for_chars(unit, chars)
    integer, intent(in) :: unit
    character(len = *), intent(in) :: chars
    integer :: i, length
    logical :: eof
    character :: c

    length = len_trim(chars)

    do i = 1, length
       c = pop_char(unit, eof=eof, skip_ws=.true.)
       if (eof) then
          print *, "ERROR: Unexpected end of file while parsing array."
          call exit (1)
       else if (c /= chars(i:i)) then
          print *, "ERROR: Unexpected character.'", c,"'", chars(i:i)
          call exit (1)
       end if
    end do

  end subroutine parse_for_chars

  !
  ! PARSE NUMBER
  !
  subroutine parse_number(unit, value)
    use fson_value_m, only: TYPE_INTEGER, TYPE_REAL
    implicit none

    integer, intent(inout) :: unit
    type(fson_value), pointer :: value
    logical :: eof, negative, decimal, scientific
    character :: c
    integer :: integral, exp, digit_count
    real(kind(1.0D0)) :: frac

    ! first character is either - or a digit
    c = pop_char(unit, eof=eof, skip_ws=.true.)
    if (eof) then
       print *, "ERROR: Unexpected end of file while parsing number."
       call exit (1)
    else if ("-" == c) then
       negative = .true.
    else
       negative = .false.
       call push_char(c)
    end if

    ! parse the integral
    integral = parse_integer(unit)

    decimal = .false.
    scientific = .false.

    do
       ! first character is either - or a digit
       c = pop_char(unit, eof=eof, skip_ws=.true.)
       if (eof) then
          print *, "ERROR: Unexpected end of file while parsing number."
          call exit (1)
       else
          select case (c)
          case (".")
             ! this is already fractional number
             if (decimal) then
                ! already found a decimal place
                print *, "ERROR: Unexpected second decimal place while parsing number."
                call exit(1)
             end if
             decimal = .true.
             frac = parse_integer(unit, digit_count)
             frac = frac / (10.0D0 ** digit_count)
          case ("e", "E")
             ! this is already an exponent number
             if (scientific) then
                ! already found a e place
                print *, "ERROR: Unexpected second exponent while parsing number."
                call exit(1)
             end if
             scientific = .true.
             ! this number has an exponent
             exp = parse_integer(unit)

          case default
             ! this is a integer
             if (decimal) then

                ! add the integral
                frac = frac + integral

                if (scientific) then
                   ! apply exponent
                   frac = frac * (10.0D0 ** exp)
                end if

                ! apply negative
                if (negative) then
                   frac = frac * (-1)
                end if

                value % value_type = TYPE_REAL
                value % value_real = frac
                value % value_double = frac
             else
                if (scientific) then
                   ! apply exponent
                   integral = integral * (10.0D0 ** exp)
                end if

                ! apply negative
                if (negative) then
                   integral = integral * (-1)
                end if

                value % value_type = TYPE_INTEGER
                value % value_integer = integral
                !+PJK
                !  Following two lines are included in case the decimal point of
                !  a floating point number has been (accidentally) left out
                value % value_real = integral
                value % value_double = integral
                !-PJK
             end if
             call push_char(c)
             exit
          end select
       end if
    end do

  end subroutine parse_number

  !
  ! PARSE INTEGER
  !
  integer(kind=8) function parse_integer(unit, digit_count) result(integral)
    integer, intent(in) :: unit
    integer, optional, intent(inout) :: digit_count
    logical :: eof, found_sign, found_digit
    character :: c
    integer :: tmp, icount, isign
    integer, parameter :: max_integer_length = 18

    icount = 0
    integral = 0
    isign = 1
    found_sign = .false.
    found_digit = .false.
    do
       c = pop_char(unit, eof=eof, skip_ws=.true.)
       if (eof) then
          print *, "ERROR: Unexpected end of file while parsing digit."
          call exit (1)
       else
          select case(c)
          case ("+")
             if (found_sign.or.found_digit) then
                print *, "ERROR: Miss formatted number."
                call exit(1)
             end if
             found_sign = .true.
          case ("-")
             if (found_sign.or.found_digit) then
                print *, "ERROR: Miss formatted number."
                call exit(1)
             end if
             found_sign = .true.
             isign = -1
          case ("0":"9")
             found_sign = .true.
             if (icount > max_integer_length) then
                print *, "ERROR: Too many digits for an integer."
                call exit(1)
             end if
             ! digit
             read (c, '(i1)') tmp
             ! shift
             if (icount > 0) then
                integral = integral * 10
             end if
             ! add
             integral = integral + tmp

             ! increase the icount
             icount = icount + 1
          case default
             if (present(digit_count)) then
                digit_count = icount
             end if
             call push_char(c)
             integral = isign * integral
             return
          end select
       end if
    end do

  end function parse_integer

  !
  ! POP CHAR
  !
  recursive character function pop_char(unit, eof, skip_ws) result(popped)
    integer, intent(in) :: unit
    logical, intent(out) :: eof
    logical, intent(in), optional :: skip_ws

    integer :: ios
    character :: c
    logical :: ignore

    !+PJK
    ios = 0
    !-PJK
    eof = .false.
    if (.not.present(skip_ws)) then
       ignore = .false.
    else
       ignore = skip_ws
    end if

    do
       if (pushed_index > 0) then
          ! there is a character pushed back on, most likely from the number parsing
          c = pushed_char(pushed_index:pushed_index)
          pushed_index = pushed_index - 1
       else
          read (unit=unit, fmt="(a)", advance="no", iostat=ios) c
       end if
       if (ios == end_of_record) then
          cycle
       else if (ios == end_of_file) then
          eof = .true.
          exit
       else if (iachar(c) <= 31) then  !  PJK from 32
          ! non printing ascii characters
          cycle
       else if (ignore .and. c == " ") then
          cycle
       else
          popped = c
          exit
       end if
    end do

  end function pop_char

  !
  ! PUSH CHAR
  !
  subroutine push_char(c)
    character, intent(inout) :: c
    pushed_index = pushed_index + 1
    pushed_char(pushed_index:pushed_index) = c

  end subroutine push_char

end module fson_library

# Profile Abstract Base Class | `Profile(ABC)`

The Profile class serves as a template for subclasses that represent different types of profiles. It initialises the profile with a specified size and sets up the basic structure for storing profile data. The profile size attribute is assigned as an integer with the profile `x` and `y` values initialised to empty numpy arrays. 

## Normalise the profile in `x` | `normalise_profile_x()`

The values of the profiles `x` dimension are normalised to be between 0 and 1 by dividing the total profile size.

## Calculate the profile steps in `x` | `calculate_profile_dx()`

The difference or step size between each concurrent value in `x` is calculated by dividing the max value in `x` by the profile size minus one.

## Calculate the profile `y` values | `calculate_profile_y()`

This acts as an abstract holder method for the particular profile solving functions to be assigned to.

## Calculate the profile integral value | `integrate_profile_y()`

The profile is integrated between its minimum and maximum bounds using the [Simpsons rule integration method](https://docs.scipy.org/doc/scipy/reference/generated/scipy.integrate.simpson.html) from `scipy` using the step size from `calculate_profile_dx()`.
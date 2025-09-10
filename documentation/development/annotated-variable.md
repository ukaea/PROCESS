# Annotating variables in Python for their inclusion in the dictionaries
As the Python conversion continues, we want to maintain backwards compatibility with our utility tools. These tools rely on the so-called 'dictionaries' to provide some of the meta data used in the utility functions.

The dictionaries, in Fortran, are created by source code analysis using `ford`. With Python, we can use the known code structure to programatically extract this information and include it within the dictionary creation - therefore maintaining backwards compatibility during the conversion work. 

## The AnnotatedVariable
This **function** can be found in `process.variables`. 

A high-level overview of this function is that it dynamically creates a class that acts like a specific type. When dynamically creating our types, we can optionally provide a `docstring` and `units`. The `docstring` provides a means of giving a variable a docstring, something not possible in Python (although most auto-documenters do allow it). The `units` are used by `AnnotatedVariable` to add a `__units__` attribute to the variable, which can be accessed by `create_dicts.py` later.

All arguments, except `tp`, and all keyword-arguments other than `docstring` and `units` are passed directly to `tp`'s constructor (`__init__`). See the example below for how this is used.

The inner working of the `AnnotatedVariable` is very abstract and, as such, will not be covered here. The important thing to understand is that the `AnnotatedVariable` function returns an object that, although claims to be of type `process.variables.AnnotatedVariable.<locals>._Variable` is actually of the type, `tp`, you provide it. 

## Using AnnotatedVariables
It only makes sense to use `AnnotatedVariable`s inside of a classes `__init__` method. This is because, the dictionary creation scripts will only look at a base initialised class, ie only the `__init__` method will be run before the object is interogated.

!!! example "Basic use of an AnnotatedVariable"

    ```python
    from process.variables import AnnotatedVariable

    class SomePhysicsModule:
        def __init__(self):
            self.var1 = AnnotatedVariable(float, 0.0)
    ```

    `self.var1` has been created as a basic `float` (standard Python type) and has been given an initial value of `0.0`. The `0.0` being an argument that is not `tp` is provided to the `float`s constructor. This code, therefore, is equivalent to the following:

    ```python
    class SomePhysicsModule:
        def __init__(self):
            self.var1 = 0.0
    ```


!!! example "More advanced use of an AnnotatedVariable"

    ```python
    from process.variables import AnnotatedVariable

    class SomePhysicsModule:
        def __init__(self):
            self.cost = AnnotatedVariable(float, 0.0, docstring="The cost associated with SomePhysicsModule", units="£")
    

    x = SomePhysicsModule()
    print(f'{type(x.cost)=}')
    print(f'{x.cost=}')
    print(f'{x.cost.__doc__=}')
    print(f'{x.cost.__units__=}')
    ```

    Here we have declated an instance variable, `self.cost` that uses both a `docstring` and `units`.

    When running this simple script, we get the following output:
    ```
    type(x.cost)=<class 'process.variables.AnnotatedVariable.<locals>._Variable'>
    x.cost=0.0
    x.cost.__doc__='The cost associated with SomePhysicsModule'
    x.cost.__units__='£'
    ```

!!! warning
    Once we act upon `self.cost`, it may lose its `__doc__` and `__units__`. That is why this idea must only be used inside of the `__init__` method.

    ```python
    from process.variables import AnnotatedVariable

    class SomePhysicsModule:
        def __init__(self):
            self.cost = AnnotatedVariable(float, 0.0, docstring="The cost associated with SomePhysicsModule", units="£")

        def do_something(self):
            self.cost += 1000.0

    x = SomePhysicsModule()
    x.do_something()
    print(f'{type(x.cost)=}')
    print(f'{x.cost=}')
    print(f'{x.cost.__doc__=}')
    print(f'{x.cost.__units__=}')
    ```

    And we get that the docstring, `__doc__` has changed and the `__units__` attribute no longer exists.

    ```
    type(x.cost)=<class 'float'>
    x.cost=1000.0
    x.cost.__doc__='Convert a string or number to a floating point number, if possible.'
    Traceback (most recent call last):
    File "temp.py", line 14, in <module>
        print(f'{x.cost.__units__=}')
    AttributeError: 'float' object has no attribute '__units__'
    ```



!!! example "Using numpy with AnnotatedVariable"

    ```python

    class SomePhysicsModule:
        def __init__(self):
            self.cost = AnnotatedVariable(np.ndarray, (5,5), docstring="The cost associated with SomePhysicsModule", units="£")

        def do_something(self):
            self.cost += 1000.0

    x = SomePhysicsModule()
    x.do_something()
    print(f'{x.cost=}')
    print(f'{x.cost.__doc__=}')
    print(f'{x.cost.__units__=}')
    ```

    We get the following output:
    ```
    x.cost=_Variable([[1000., 1000., 1000., 1000., 1000.],
           [1000., 1000., 1000., 1000., 1000.],
           [1000., 1000., 1000., 1000., 1000.],
           [1000., 1000., 1000., 1000., 1000.],
           [1000., 1000., 1000., 1000., 1000.]])

    x.cost.__doc__='The cost associated with SomePhysicsModule'

    x.cost.__units__='£'
    ```

!!! note
    Note in the above example that `x.cost` is shown as type _Variable. This is to do with the way the AnnotatedVariable wrapper injects data into a class.

    However, checking that `x.cost` is a numpy array confirms that this is indeed still a numpy array:

    ```python
    print(f'{isinstance(x.cost, np.ndarray)=}')
    ```

    ```
    isinstance(x.cost, np.ndarray)=True
    ```

!!! note
    In this instance, despite operating on the data structure in the `do_something` method, we have retained the `__doc__` and `__unit__` methods.

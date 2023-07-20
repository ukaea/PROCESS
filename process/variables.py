"""Provides wrappers around variables to the whole of PROCESS that hydrate the code with various functionalities"""

from typing import Type


def AnnotatedVariable(
    tp: Type, *args, docstring: str = "", units: str = "", __kwargs: dict = {}, **kwargs
):
    """Provides a wrapper around the instantiation of a variable to allow variables to be type hinted. This should be done on physics and engineering class variables to allow for their inclusion in 'the Dictionaries'.

    The returned object is a subclass of `tp` and can be operated on as such. The only difference, is that this object has `__doc__` and `__units__` attributes, corresponding, respectively, to `docstring` and `units`.

    :param tp: the type of variable to create e.g. float, int
    :type tp: Type

    :param *args: any positional arguments to be passed to the constructor of the type

    :param docstring: the docstring of the variable
    :type docstring: str

    :param units: the units of the variable e.g. kg
    :type units: str

    :param __kwargs: provides a means to specify kwargs to the constructor of the type that clash with names in this functions parameters e.g. tp(units='') __kwargs={'units': 'same keywork, different value and purpose'}
    :type __kwargs: Dict[str, Any]

    :param **kwargs: keyword arguments provided to the constructor of the type
    """
    if not isinstance(tp, type):
        raise TypeError(f'Argument type_ must be of type "type" not {type(tp)}')

    class _Variable(tp):
        __doc__ = docstring
        __units__ = units

    return _Variable(*args, **{**kwargs, **__kwargs})

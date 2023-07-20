from django import template

register = template.Library()


@register.filter
def get_item(dictionary, key):
    """Gets an item from the dictionary.
    rather than [] so will return None when key
    does not exist
    """
    return dictionary[key]


@register.filter
def to_str(x):
    """Converts an item to string format
    Lists are stripped of surrounding [ ] characters
    Strings are given surrounding ' ' quotes
    """
    if isinstance(x, list):
        return str(x)[1:-1]
    elif isinstance(x, str):
        return "'" + str(x) + "'"
    else:
        return x


@register.filter
def tolbname(x):
    """Converts a ixc no to  boundl(#)"""
    return "boundl(" + str(x) + ")"


@register.filter
def toubname(x):
    """Converts a ixc no to  boundu(#)"""
    return "boundu(" + str(x) + ")"

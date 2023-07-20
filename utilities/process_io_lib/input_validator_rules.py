"""The rules for the input validator to run.

This module contains the rule classes which are used by the input_validator
module to check rules on the input data to validate it. This consists of the
abstract Rule class and the individual rule subclasses which are intialised by
the input_validator module. Each rule class (a subclass of Rule) defines its own
check method, which checks that particular rule against the input data. This
method stores the result and any messages on the instance of that rule class
itself.

To add a new rule, define a new class that inherits from the Rule class, naming
it the same as the variable it covers, but with the first letter capitalised.
Then override the __init__() and check() methods. The new rule class will be
used by the input_validator module automatically, and the rule will be checked
when input_validator is run. See class Ishape(Rule) for an example, or use the
rule snippet in the Process project on Gitlab.
"""
from abc import ABC, abstractmethod
from process_io_lib.obsolete_vars import OBS_VARS


class Rule(ABC):
    """Abstract rule class used by individual rule subclasses

    Each rule to check on the input file data is a subclass of this. This is an
    abstract base class, so only subclasses of this can be instantiated.
    """

    @abstractmethod
    def __init__(self, tags):
        """Initialise a Rule object.

        This is an abstract method: it must be overridden in subclasses before
        they can be instantiated.

        :param tags: Tags for describing the rule; used for filtering rules
        (list of strings)
        :type tags: list
        """
        self.name = self.__class__.__name__.lower()
        # Set the rule's name to the name of the class. Lower to match the
        # variable name that the rule checks
        self.tags = tags
        self.passed = True
        # Stores result of rule. Set to False in check method in case of
        # failure
        self.messages = []
        # List of strings storing reason(s) for check failure
        self.data = None
        # Data for the rule to run on

    def set_data(self, data):
        """Save the input data onto the Rule object so the check method has
        access to it.

        :param data: The input data object
        :type data: object
        """
        self.data = data

    @abstractmethod
    def check(self):
        """Method to check the rule, and set the passed attribute.

        Decides whether the passed attribute should be set to True or False.
        This is an abstract method: it must be overridden in subclasses before
        they can be instantiated.
        """
        pass

    def check_defined(self, var_name):
        """Checks that a parameter is defined.

        If the parameter is undefined, sets the self.passed attribute to False
        and stores a message to record the failure of the check.

        :param var_name: The name of the parameter
        :type var_name: str
        """
        defined = self.data.is_param_defined(var_name)

        if not defined:
            message = f"{var_name} should be defined."
            self.messages.append(message)
            self.passed = False

    def check_undefined(self, var_name):
        """Checks that a parameter is undefined.

        If the parameter is defined, sets the self.passed attribute to False
        and stores a message to record the failure of the check.

        :param var_name: The name of the parameter
        :type var_name: str
        """
        defined = self.data.is_param_defined(var_name)

        if defined:
            message = f"{var_name} shouldn't be defined."
            self.messages.append(message)
            self.passed = False

    def get_param_value(self, var_name):
        """Gets the value of a parameter in the input data.

        :param var_name: Name of the parameter
        :type var_name: str
        :return: Value of the parameter
        :rtype: int, float
        """
        return self.data.get_param_value(var_name)


# Rule classes for checking individual rules
class Ishape(Rule):
    """Rule subclass for checking the value of ishape and its dependencies

    ishape: switch for plasma cross-sectional shape calculation.
    """

    def __init__(self):
        """Call Rule's __init__ method with tags specific to Ishape"""
        super().__init__(["ishape", "kappa", "triang", "kappa95", "triang95"])

    def check(self):
        """The rule function for Ishape to check on the input data"""
        # ishape must be defined
        self.check_defined("ishape")
        ishape = self.get_param_value("ishape")

        if ishape == 0:
            # Use kappa and triang to calculate 95% kappa and triang
            self.check_defined("kappa")
            self.check_defined("triang")
            self.check_undefined("kappa95")
            self.check_undefined("triang95")
        elif ishape == 1:
            # Scale with aspect ratio
            self.check_undefined("qlim")
            self.check_undefined("kappa")
            self.check_undefined("triang")
            self.check_undefined("kappa95")
            self.check_undefined("triang95")
        elif ishape == 2:
            # kappa calculated using fkzohm, triang input
            self.check_defined("fkzohm")
            self.check_defined("aspect")
            self.check_undefined("kappa")
            self.check_defined("triang")
            self.check_undefined("kappa95")
            self.check_undefined("triang95")
        elif ishape == 3:
            # kappa calculated using fkzohm, triang95 input
            self.check_defined("fkzohm")
            self.check_undefined("kappa")
            self.check_defined("triang95")
            self.check_undefined("triang")
            self.check_undefined("kappa95")
        elif ishape == 4:
            # kappa95 and triang95 are used to calculate kappa and triang
            self.check_defined("kappa95")
            self.check_defined("triang95")
            self.check_undefined("kappa")
            self.check_undefined("triang")
        elif ishape == 5:
            # kappa95 and triang95 are used to calculate kappa and triang
            self.check_defined("kappa95")
            self.check_defined("triang95")
            self.check_undefined("kappa")
            self.check_undefined("triang")
        elif ishape == 5:
            # Use kappa and triang to calculate 95% kappa and triang
            self.check_defined("kappa")
            self.check_defined("triang")
            self.check_undefined("kappa95")
            self.check_undefined("triang95")
        elif ishape == 7:
            # kappa95 and triang95 are used to calculate kappa and triang
            self.check_defined("kappa95")
            self.check_defined("triang95")
            self.check_undefined("kappa")
            self.check_undefined("triang")
        elif ishape == 8:
            # Use kappa and triang to calculate 95% kappa and triang
            self.check_defined("kappa")
            self.check_defined("triang")
            self.check_undefined("kappa95")
            self.check_undefined("triang95")


class Aspect(Rule):
    """Aspect ratio"""

    def __init__(self):
        """Set tags specific to Aspect"""
        super().__init__(["aspect"])

    def check(self):
        """Check that aspect exists in input"""
        self.check_defined("aspect")


class Hfact(Rule):
    """Energy confinement time H-factor"""

    def __init__(self):
        """Set tags specific to Hfact"""
        super().__init__(["hfact"])

    def check(self):
        """Check hfact is defined"""
        self.check_defined("hfact")


class Kappa(Rule):
    """Plasma elongation"""

    def __init__(self):
        """Set tags specific to Kappa"""
        super().__init__(["kappa", "ishape"])

    def check(self):
        """Should kappa be defined, based on value of ishape"""
        # This logic is covered in the Ishape rule
        ishape = self.get_param_value("ishape")
        if ishape is [0, 6, 8]:
            # kappa input value is used
            self.check_defined("kappa")
        else:
            # kappa is calculated
            self.check_undefined("kappa")


class Triang(Rule):
    """Plasma triangularity"""

    def __init__(self):
        """Set tags specific to Triang"""
        super().__init__(["triang", "ishape"])

    def check(self):
        """Should triang be defined, based on value of ishape"""
        # This logic is covered in the the Ishape rule
        ishape = self.get_param_value("ishape")
        if ishape in [0, 2, 6, 8]:
            self.check_defined("triang")
        else:
            self.check_undefined("triang")


class Alphan(Rule):
    """Density profile index"""

    def __init__(self):
        """Set tags specific to Alphan"""
        super().__init__(["alphan"])

    def check(self):
        """Check alphan is defined"""
        self.check_defined("alphan")


class Alphat(Rule):
    """Temperature profile index"""

    def __init__(self):
        """Set tags specific to Alphat"""
        super().__init__(["alphat"])

    def check(self):
        """Check alphat is defined"""
        self.check_defined("alphat")


class Dnbeta(Rule):
    """(Troyon-like) coefficient for beta scaling"""

    def __init__(self):
        """Set tags specific to Dnbeta"""
        super().__init__(["dnbeta", "gtscale"])

    def check(self):
        """Check if dnbeta input value is required or not"""
        iprofile = self.get_param_value("iprofile")
        gtscale = self.get_param_value("gtscale")
        if iprofile == 0:
            if gtscale == 1:
                # dnbeta is calculated
                self.check_undefined("dnbeta")
            else:
                # dnbeta is required
                self.check_defined("dnbeta")
        else:
            # dnbeta is calculated
            self.check_undefined("dnbeta")


class ObsoleteVarChecker(Rule):
    """Checks for obsolete or renamed vars, and suggests new names."""

    def __init__(self):
        """Set tags specific to ObsoleteVarChecker"""
        super().__init__(["obsolete"])

    def check(self):
        """Look for obsolete variable names and suggest alternatives."""
        # List of unrecognised vars in the input file
        unrecog_vars = self.data.unrecognised_vars

        # If there are obsolete vars in the input file, fail the rule, warn
        # and suggest alternatives, if available
        for unrecog_var in unrecog_vars:
            # Try to find the unrecognised var in the obsolete vars dict
            recog_var = OBS_VARS.get(unrecog_var, False)

            if recog_var:
                # Obsolete var name with new var name suggestion found
                message = f'"{unrecog_var}" is obsolete; try "{recog_var}"' " instead."
            elif recog_var is None:
                # Obsolete var name, no new var name suggestion
                message = f"{unrecog_var} is deprecated."
            else:
                # Var name is unrecognised, but not in obsolete list either
                break

            # Fail rule with message
            self.passed = False
            self.messages.append(message)


class DuplicateChecker(Rule):
    """Ensures there are no duplicate variable initialisations."""

    def __init__(self):
        """Set tags specific to DuplicateChecker"""
        super().__init__(["duplicates"])

    def check(self):
        """Find any duplicate variables in the input data."""
        duplicates = self.data.duplicates
        if duplicates:
            # Got some duplicates
            self.passed = False
            duplicates_no = len(duplicates)
            message = f"Found {duplicates_no} duplicates: {duplicates}"
            self.messages.append(message)

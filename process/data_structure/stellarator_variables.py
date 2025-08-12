# These variables were from stellarator.f90
f_n: float = None

f_r: float = None

f_aspect: float = None

f_b: float = None

f_i: float = None

f_a: float = None

first_call: bool = None

first_call_stfwbs: bool = None


def init_stellarator_module():
    global first_call
    global first_call_stfwbs
    global f_n
    global f_r
    global f_a
    global f_b
    global f_i

    first_call = True
    first_call_stfwbs = True
    f_n = 0.0
    f_r = 0.0
    f_a = 0.0
    f_b = 0.0
    f_i = 0.0

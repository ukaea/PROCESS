import os
import re


def system_compatible():
    """Check the system compatibility to ensure tests are not
    overwritten with inaccurate data caused by floating point differences in
    dynamically linked libraries provided by the OS.

    :return: is the system compatible with the test suite, and can it update the test
    assets?
    :rtype: boolean
    """
    try:
        shell_stream = os.popen("ldd --version")

        ldd_version = re.search(r"^ldd .+ (2.[0-9.]+).*", shell_stream.read()).group(1)
        ldd_version_array = ldd_version.split(".")
        major = ldd_version_array[0]
        minor = ldd_version_array[1]

        if int(major) < 2 or (int(major) == 2 and int(minor) < 31):
            return False
    except AttributeError:
        return False

    return True

# Introduction

Below are instructions of how to install the applications, code and dependencies such that you
can begin working on PROCESS. This includes differing instructions for those working on
different operating systems.

!!! Info "Python version"
    Please note, only Python3.8 is supported, Python3.9, Python3.10 and Python3.6 are not.

PROCESS is natively supported on Ubuntu 20. Other Linux distributions will be able to successfully
build and execute PROCESS however may give inconsistent results due to version disparities of
dynamically linked libraries.

There are three supported ways to run PROCESS on a non-native Ubuntu 20.04 machine:

1. WSL (Windows Subsystem for Linux- if you are a Windows user) --> [guide](installation-ubuntu.md)
2. Docker (for users on privileged machines that are of the wrong OS) --> [guide](installation-docker.md)
3. Singularity (mainly for use on shared resources e.g. Freia) --> [guide](installation-singularity.md)

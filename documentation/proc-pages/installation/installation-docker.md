# Installation using Docker

!!! Note "Note"
    This section is for users of a privileged machine of the wrong OS e.g MacOS. If you are a using a shared resource e.g. computing cluster, please go [here](#singularity-container).*

PROCESS can be run on Mac or in other environments inside a Docker container. The PROCESS 
repository, including source and build directories, remain in the host filesystem, but the 
building and running of PROCESS is performed inside the container. This ensures that PROCESS 
produces the same results as in other fully-supported environments, such as the CI system. The 
Ubuntu-based development image used is similar to the one used on the CI system, but it is 
designed to work immediately with no further installations.

!!! Note "Docker licence"
    Please note due to recent changes in the Docker Desktop ToS, you will require either a Docker 
    Desktop license to run on Mac, or you will require a Linux environment by other means, such 
    as a virtual machine.

Firstly, [install Docker](https://docs.docker.com/get-docker/). On Mac, this can be accomplished 
using `homebrew`:

```
brew --cask install docker
```

Then login to the Gitlab container registry:

```
docker login git.ccfe.ac.uk:4567
```

Then download the Docker image from the Process Gitlab container registry:

```
docker pull git.ccfe.ac.uk:4567/process/process/dev
```

Running `docker image ls` should show the image in your local Docker image repository. Optionally, 
you can change the image name to something more manageable:

```
docker tag git.ccfe.ac.uk:4567/process/process/dev process-dev
```

to rename the image to "process-dev" with the "latest" tag: "process:latest".

Now run the container:

```
docker run -it -v ~/process:/root/process process-dev
```

This runs a container which is an instance of the process-dev image. `-it` runs the container in 
interactive mode (`-i`, allows `stdin`) with a terminal (`-t`, allows bash-like interaction). `-v` 
specifies the bind mount to use; mount the host `~/process` directory to the `/root/process` 
directory in the container. This means that the container has read and write access to the `process` 
project directory on the host filesystem and will stay in sync with it. Please be aware that 
changes made in a Docker container outside of mounted folders will not be saved on exiting the container.

Now the container is running, configure, clean and build from the project root directory inside the container:

```
cd ~/process
cmake -S . -B build
cmake --build build --clean-first
```

The clean step is required to remove any build targets or caches from previous host builds to ensure 
a completely fresh build from inside the container. This is only required when using the container 
for the first time.

Once PROCESS has built inside the container, it can be tested (as in the following section) by 
running `pytest`. Once the test suite passes, this confirms that your Docker container runs PROCESS 
with the same results as the CI system. PROCESS can now be developed and run as before, with the 
build and running taking place inside the container.

There is also a VS Code extension for Docker containers that may be helpful.

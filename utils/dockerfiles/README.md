BubbleProfiler Dockerfiles
==========================

The following Dockerfiles are provided for
testing BubbleProfiler:

  * bubbleprofiler-base-ubuntu-gcc
    - contains a base Ubuntu image with gcc and
      make installed
  * bubbleprofiler-ci-ubuntu-gcc
    - based on bubbleprofiler-base-ubuntu-gcc,
      this image also installs the necessary
      dependencies for building and running
      BubbleProfiler
  * bubbleprofiler-ci-ubuntu-docs
    - based on bubbleprofiler-ci-ubuntu-gcc,
      this image also installs Doxygen for
      the purpose of generating the API
      documentation

Update procedure
----------------

To update an image in the public Docker
repository,

1. First build the image locally with
  ```bash
  cd <image directory>
  docker build -t <image name> .
  ```

2. Log-in to the public Docker repository,
  ```bash
  docker login
  ```

3. Tag the update image appropriately,
  ```bash
  docker tag <image name> bubbleprofiler/<image name>:devel
  docker tag <image name> bubbleprofiler/<image name>:latest
  ```

4. Push the updated image to the public Docker repository,
  ```bash
  docker push bubbleprofiler/<image name>:devel
  docker push bubbleprofiler/<image name>:latest
  ```

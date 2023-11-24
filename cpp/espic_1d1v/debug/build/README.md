# Summary
`Eigen` must be installed, and CMake must be pointed to the folder where the `FindEigen3.cmake` file is. 

On my machine, this is accomplished via,

    cmake ../src -DCMAKE_MODULE_PATH=$HOME/software/eigen-3.4.0/cmake

Generally, it will be of the form,

    cmake ../src -DCMAKE_MODULE_PATH=$HOME/software/eigen-3.4.0/cmake

# toptor

Introduction
============
The TopTor tool has been designed to call the TopGen Library from the  command line interface. It access to all the functionality of the library without the necessity of having a simulator. 

Get Started
===========

How to compile
--------------

The first step to compile the tool consists in cloning the repository and invoking cmake with the appropriate parameters.
 
```bash

DIR="$HOME/opt"

cmake -DCMAKE_INSTALL_PREFIX="$DIR/toptor" -DTOPGEN="$DIR/topgen" -DBOOST_ROOT="$DIR/boost-165" -DCMAKE_BUILD_TYPE=Debug .
make -j 5
make install
```
If the compilation is successful, then the tool toptor will be inside `$DIR/toptor/bin`.

Credits
=======

The authors

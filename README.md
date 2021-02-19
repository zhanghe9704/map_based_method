## Map_based_method lib

### Announcement

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

### About the code

This code includes some map-based methods for accelerator and beam dynamic simulation and analysis. However, at least for now, this code does NOT generate the map. The users have to define their system and generate the map using COSY Infinity or MAD-X. This code can read the output file from them. The [tpsa lib](https://github.com/zhanghe9704/tpsa) is also included. No need to load the tpsa lib separately for Differential Algebra (DA) or  Truncated Power Series Algebra (TPSA) calculations. 

By 02/19/2021, the following tools are included:

Tracking using a truncated Taylor map

Symplectic tracking using the generating functions derived from a truncated Taylor map [1]

Lie factorization for a truncated Taylor map [2]



> [1] **Symplectic tracking in circular accelerators with high order maps** by Martin Berz in "Nonlinear Problems in Future Particle Accelerators" (1991) 288-296, World Scientific. 
>
> [2] **Arbitrary order description of arbitrary particle optical systems** by Martin Berz in "Nuclear Instruments and Methods in Physics Research A298" (1990) 426-440



### How to compile

The following process is tested in Ubuntu 18.04 in Windows 10 Subsystem for Linux (WSL) environment. 

To compile the source file, a gcc compiler that supports C++ 11 is needed. A Makefile is provided. Assume the repository is cloned into $HOME/map_based_method. Enter this directory and run 

```shell
make
```

A shared lib named "libmap_method.so" will be created. 

### How to compile and run the examples

Enter the directory $HOME/map_based_method/examples and run

```shell
make
```

Executable files for all the examples will be created.  To run the examples, the shared lib "libmap_method.so" is needed, which is created by the source code (see above).  Copy the lib into the "examples" folder. To make sure they can be found by the system, in the "examples" folder, run

```shell
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$PWD
```

Mow you ready to run the executable files. For example, with the lib and the input file ready, one can run example_tracking_01 by  

```shell
./example_tracking_01
```



The example\_tracking_01 reads a 4D map generated by COSY Infinity 9.x and then applies both the map tracking and symplectic tracking. "4D\_2ndOrder_COSY.txt" is the input file for example\_tracking\_01. "track\_truncated\_map\_cosy.txt" and "track\_symplectic\_type\_1\_cosy.txt" are the outputs from example\_tracking\_01. 

The example\_tracking_02 reads a 6D map generated by MAD-X and then applies both the map tracking and symplectic tracking. "6D\_2ndOrder_MADX.txt" is the input file for example\_tracking\_02. "track\_truncated\_map\_madx.txt" and "track\_symplectic\_type\_2\_madx.txt" are the outputs from example\_tracking\_02. 

The file "track.plt" plots the results using gnuplot. 

The example\_tracking_03 shows how to separate the symplectic tracking process into two steps: (1) derive the differential equation that enforces the symplecticity and (2) perform the syplectic tracking by solving the differential equation for each particle. 

The example\_tracking_04 shows how to convert a map that takes PTC coordinates into one that takes MAD-X coordinates and use the maps for tracking. 

The example_lie_factorization shows how to perform Lie factorization on a truncated Taylor map. 

### Contact

Contact the author by hezhang.AT.jlab.org.


# GraphCut
A C++ implement of the paper "Graphcut Textures: Image and Video Synthesis Using Graph Cuts"

## Compile
```bash
mkdir build
cd build
cmake .. -DCMAKE_BUILD_TYPE=Release
make
```

## Usage

```
# Example: graph_cut peas.png peas_output.png 512x512
graph_cut <input> <output> <canvas_size>
```

## Details

For more details, please refer to report.pdf.



## Performance support setup

Need to install OpenMP and FFTW on the system.
Using FFTW-3.3.9

```
#To install fftw:
tar -xvf fftw-3.3.9.tar	
cd fftw-3.3.9.tar
./configure
make
make install
```

.


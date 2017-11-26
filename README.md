# morph
This is my personal research code for morphing multiple genus zero shapes. Code is written in Matlab. Use at your own risk. 

## How to use the code
My codes are under `morphing_codes_max`. Put the mesh files under `mesh_files`. Run and go through `run_demo.m` for an example on how to use the code. The main code is `morphing_codes_max/morph_multi.m`, which returns a face matrix `F` and various vertex matrix `VS1`, `VS2`, ... Create a morph by linearly interpolating the vertex matrices. Then you may export it to a standard mesh format. e.g
```matlab
n_feat = 9; % Number of feature points to match between models
filename1 = 'cow40k.ply';
filename2 = 'horse50k.ply';
[F,VS1,VS2] = morph_multi(n_feat,filename1,filename2);
 V = 0.5 * VS1 + 0.5 * VS2;
 writeSTL('morphed_mesh.stl', V, F);
 ```

## Compile mex code in gptoolbox
In general, refer to [alecjacobson/gptoolbox](https://github.com/alecjacobson/gptoolbox) for instruction on compiling gptoolbox (which is included in a subdirectory in morph). The compilation guidelines apply to Mac/Linux computers. For Windows computers, refer to the getaround below.

## Windows Users
### Compile mex code
1. In Matlab go into the folder `gptoolbox/mex`
Then enter the following:
```matlab
MEXOPTS={'-v','-largeArrayDims','-DMEX'};
STDCPP11='CXXFLAGS=$CXXFLAGS -std=c++11';
mex( ...
MEXOPTS{:}, ...
STDCPP11, ...
'ray_mesh_intersect.cpp');
```

2. If that doesn't work, it's probably because your computer does not have a SDK that contains a C compiler. In MATLAB, go to `Home -> Environment -> Add-ons ->` Search for `"MinGW64"` and then install the add-on. Then try (1) again.

3. Move the generated `ray_mesh_intersect.mexw64` to `gptoolbox/mesh`
### Set path
In the demo code `run_demo` line 4, change `system('pwd')` to `system('cd')` in Windows.

# morph
Code to morph multiple genus zero shapes (Matlab)

## Compile mex code in gptoolbox
In general, refer to [alecjacobson/gptoolbox](https://github.com/alecjacobson/gptoolbox) for instruction on compiling gptoolbox (which is included in a subdirectory in morph). The compilation guidelines apply to Mac/Linux computers. For Windows computers, refer the the getaround below.
### Windows
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

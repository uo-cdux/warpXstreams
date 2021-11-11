# Prerequisite

For generating WarpX streamlines using this code you'll require:
1. Boost (program options)
2. VTK-m version 1.7

Sample WarpX data is available [here](https://www.dropbox.com/s/z2m4psgjqetxqgl/warpXdata.tar?dl=0)
There are two different types of data:
1. Species : These are the charged particles that will be used for advection.
2. Fields  : These are the electric and magnetic fields that we'll use for velocity calculation.

To execute the code after compilation you can simply execute the following:
```
./advection params
```
`params` being the file containing the input parameters for the example.
The parameters that are included for the data at the link above are
```
data=data/vtk_fields_0000250.vtk                                                
steps=120                                                                       
length=0.1                                                                      
seeds=50                                                                        
seeddata=data/vtk_specie_beam_0000250.vtk                                       
sampleZ=-6.5000e-05:-5.00668e-05                                                
```

# Warp X data

The data in the section above is only a single slice,
More WarpX data is available [here](https://www.dropbox.com/s/nfx3z35d916miw5/2020_11_15_rotating_beam-20201117T025553Z-001.zip?dl=0)

All the data is available as `hdf5` files.
They need to be converted into `vtk` in order to use the code example from this repository.

The conversion can be done as:

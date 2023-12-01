# dielectricInterfacesGeneralized
This code implements the algorithm for applying the method of images to an arbitrary number of conducting/dielectric slabs described in [Haji-Akbari A, Shoemaker B. Ideal conductor/dielectric model (ICDM): A generalized technique to correct for finite-size effects in molecular simulations of hindered ion transport. ChemRxiv. Cambridge: Cambridge Open Engage; 2023 doi:10.26434/chemrxiv-2023-1mgbh]

The basic usage is demonstrated in the function "testFunctionality.py" and described here.

__d = buildDomain([1000000,1.2,1.3],[0,2],[1],[1],10)_
This constructs a 3-region system with a single charge in the middle layer. 
[1000000,1.2,1.3] gives the dielectric constant in each region
[0,2] gives the location of the interfaces between regions
[1] specifices a point charge added at position z=1 (this list can contain multiple values to add multiple charges)
[1] specifies the point has a charge=1 (this list can contain multiple values to add multiple charges)
10 specifies that the algorithm will undergo 10 iterations in which previous image charges undergo reflection

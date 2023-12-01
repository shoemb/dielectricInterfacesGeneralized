# dielectricInterfacesGeneralized
This code implements the algorithm for applying the method of images to an arbitrary number of conducting/dielectric slabs described in [Haji-Akbari A, Shoemaker B. Ideal conductor/dielectric model (ICDM): A generalized technique to correct for finite-size effects in molecular simulations of hindered ion transport. ChemRxiv. Cambridge: Cambridge Open Engage; 2023 doi:10.26434/chemrxiv-2023-1mgbh]

The basic usage is demonstrated in the function "testFunctionality()" and described here.

--------------------------------------------------------------

**_d = buildDomain([1000000,1.2,1.3],[0,2],[1],[1],10)_**

This constructs a 3-region system with a single charge in the middle layer. 
[1000000,1.2,1.3] gives the dielectric constant in each region.
[0,2] gives the location of the interfaces between regions.
[1] specifices a point charge added at position z=1 (this list can contain multiple values to add multiple charges).
[1] specifies the point has a charge=1 (this list can contain multiple values to add multiple charges).
10 specifies that the algorithm will undergo 10 iterations in which previous image charges undergo reflection.

-------------------------------------------------------------

**_d.checkSolutionAccuracy(.00001)_**

This numerically evaluates boundary conditions at interfaces and prints the results. The number of iterations can be increased until the conditions at each side of an interface reach acceptable convergence.

--------------------------------------------------------------

**_[q,x,y,z] = d.returnImageCharges()_**

**_plt.plot(z[1],q[1],'bo')_**

**_for region in d.regions:_**

**_&emsp;if(not isinf(region.zleft)):_**
    
**_&emsp;&emsp;plt.axvline(x=region.zleft,linestyle=':')_**
        
**_plt.xlabel('Position')_**

**_plt.ylabel('Charge')_**

**_plt.show()_**

Retrieve and plot the position and charge of all image charges and free charges which determine the potential in region 1 (numbering begins with 0, so this is the central region).

--------------------------------------------------------------


**_[z,v] = d.getPotentialFull(30)_**

**_plt.plot(z,v,'b')_**

**_plt.xlabel('Position')_**

**_plt.ylabel('Potential')_**

**_plt.show()_**

Calculate and plot the potential of the full system.

--------------------------------------------------------------


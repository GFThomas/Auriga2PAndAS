# Auriga2PAndAS

This GitHub is dedicated to the Auriga2PAndAS pipeline, described in Thomas G. F. et al 2020 (xxxxx)

## Last update: 26 January 2021 by Guillaume F. Thomas *guillaume.thomas.astro .at. gmail.com*


-------------------------------------------------------------------------------------------------------------------------------
### Access to the post-processed mocks:

For each of the 6 halos treated in Thomas et al. 2020 (Au 6, 16, 21, 23, 24, 27) , a PAndAS like mock is disponible in:  *https://wwwmpa.mpa-garching.mpg.de/auriga/gaiamock.html*.
The column of these catalogues are disbribed in the Annexe 1 section of Thomas G. F. et al 2020 (xxxxx)



-------------------------------------------------------------------------------------------------------------------------------
### Create new mocks:

It is possible to create new mocks using the *exec.py* program. Below are listed the requirements and the procedure to run the Auriga2PAndAS program:


_______________________________________________________________
#### Requirements:

* `python` 3.7 with the following packages:
            - `h5py`     (https://www.h5py.org/)
            - `numpy`  (https://numpy.org/)
            - `pandas` (https://pandas.pydata.org/)
            - `astropy` (https://www.astropy.org/)
            - `ctypes`   (https://docs.python.org/3/library/ctypes.html)

* `gcc` GNU compiler 

_______________________________________________________________
#### Installation:

	>make

This will compile the parts of the Auriga2PAndAS pipeline written in c that generate the foreground contamination and apply the PAndAS footprint. It will also download the original catalogues, so it might take a while.


_______________________________________________________________
#### Execution:

	>python exec.py catalogues/Auriga_L3_halo_xx_stars.hdf5

This will execute the Auriga2PAndAS pipeline to the Halo xx of Auriga, stored in the repository catalogues. 6 halos are available: Au 6, 16, 21, 23, 24, 27 (level resolution 3)

The final file will be :
*./results/mock_L3_halo_xx.fits*

*Attention:* this will erase the previous file, therefore, if you want to keep it, we strongly suggest that you change the name of the previous mock before launching the  Auriga2PAndAS pipeline.

_______________________________________________________________
#### Options:

The code currently a random point of view around the simulted galaxy disk, and will then project it on the sky. However, a specific angle can be specify by specifying a value line 597 (angle are in radians)

	>theta=random.uniform(0, 2.0*np.pi)


-------------------------------------------------------------------------------------------------------------------------------
### Other codes:

Their is 2 others code available in the pipeline,  *fill_gaps.py* and *extract_param.py*.

_______________________________________________________________
#### fill_gaps.py

This code will fill the gap of the PAndAS footprint by duplicating neighboor regions in the mock.
** Execution: 

	>python fill_gaps.py results/mock_L3_halo_xx.fits 
	
** Results store in *./results/mock_L3_halo_xx_gaps.fits*


_______________________________________________________________
#### extract_param.py

This code will extract  one paramater of the simulations not keep by the Auriga2PAndAS pipeline (like the age of a stars). 
The code is currently extracting the age, it can be change in lines 208 and 209.

** Execution: 

	>python extract_param.py catalogues/Auriga_L3_halo_xx_stars.hdf5 results/accreted/mock_L3_halo_xx.fits
	
** Results store in *./results/mockwparam.fits*

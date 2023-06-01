## 3D heterogeneous distribution
Reconstruction of an heterogeneous distribution of particles generated with a Random Walk Particle Tracking solver on a three dimensional heterogeneous aquifer.
 
Run reconstruction
```
gpkde gpkde.sim
```

The program writes density to `gpkde.out`. The `python` routine `exportovtk.py` generates the file `gpkde.vts` for visualization with `paraview`. The exported file contains both the histogram and smoothed density reconstruction. Be sure that the packages `numpy, pandas and pyevtk` are installed. If not, a `requirements.txt` file is included and can be installed with

```
pip install -r requirements.txt
```

Note: while modifying the grid input parameters in `gpkde.sim` like the bin or domain sizes, be sure to update these values also on `exportovtk.py` for consistency while export results.

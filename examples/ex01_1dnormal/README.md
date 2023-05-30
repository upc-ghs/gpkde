## 1D Gaussian distribution 
Source particles are selected stochastically from a standard Gaussian distribution. 

The `python` routine `plotoutput.py` generates a figure comparing the smoothed and histogram density reconstruction versus the analytical solution. Be sure that the packages `numpy, scipy, pandas and matplotlib` are installed. If not, a `requirements.txt` file is included and can be installed with

```
pip install -r requirements.txt
```

Note: while modifying the grid input parameters in `gpkde.sim` like the bin or domain sizes, be sure to update these values also on `plotoutput.py` for consistency while plotting results.

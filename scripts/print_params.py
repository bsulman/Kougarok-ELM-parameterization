import sys
import xarray
from kougarok_plotting import pft_params

paramfile=sys.argv[1]
params=sys.argv[2:]

print(pft_params(xarray.open_dataset(paramfile),params))


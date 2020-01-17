import xarray
import sys
from kougarok_plotting import pft_params

def comp_params(params1,params2):
    pdiff=params2-params1
    vars_diff=[]
    for var in pdiff.variables:
        if pdiff[var].dtype.kind in ['S','m']:
            continue
        if pdiff[var].dims==('pft',):
            if pdiff[var].dropna(dim='pft').any():
                vars_diff.append(var)
    return (pft_params(params1,vars_diff),pft_params(params2,vars_diff),pft_params(pdiff,vars_diff))


def comp_params_xarray(params1,params2,tol=1e-7):
    pdiff=(params2-params1)
    psame=(params2==params1).all()
    drop_vars=[]
    for var in psame.variables:
        if psame[var] or pdiff[var].dtype.kind in ['S','m']:
            drop_vars.append(var)
        elif 'pft' in pdiff[var].dims and not pdiff[var].dropna(dim='pft').any():
            drop_vars.append(var)
        elif (abs(pdiff[var])/((params1[var]+params2[var])/2)).to_masked_array().max()<tol:
            drop_vars.append(var)
    return params1.drop_vars(drop_vars),params2.drop_vars(drop_vars),pdiff.drop_vars(drop_vars)


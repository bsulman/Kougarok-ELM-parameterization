from pylab import *
import xarray
import numpy
import pandas

params_data=xarray.open_dataset('param_files/clm_params_newpfts_c180524_orig.nc',autoclose=True)
params_default_data=xarray.open_dataset('param_files/clm_params_defaultpfts_c180524_orig.nc',autoclose=True)

pft_names=[name.strip() for name in params_data['pftname'].values.astype(str)]
pft_names_default=[name.strip() for name in params_default_data['pftname'].values.astype(str)][:17]

PFT_mappings={'not_vegetated':'not_vegetated',
 'arctic_lichen':'c3_arctic_grass',
 'arctic_bryophyte':'c3_arctic_grass',
 'arctic_evergreen_shrub_dwarf':'broadleaf_evergreen_shrub',
 'arctic_evergreen_shrub_tall':'broadleaf_evergreen_shrub',
 'arctic_deciduous_shrub_dwarf':'broadleaf_deciduous_boreal_shrub',
 'arctic_deciduous_shrub_low':'broadleaf_deciduous_boreal_shrub',
 'arctic_deciduous_shrub_tall':'broadleaf_deciduous_boreal_shrub',
 'arctic_deciduous_shrub_alder':'broadleaf_deciduous_boreal_shrub',
 'arctic_forb':'c3_arctic_grass',
 'arctic_dry_graminoid':'c3_arctic_grass',
 'arctic_wet_graminoid':'c3_arctic_grass'}

for pftname in pft_names:
    print('\n\nArctic PFT: {arcticname:s}. Default PFT: {defaultname:s}'.format(arcticname=pftname,defaultname=PFT_mappings[pftname]))
    diff_params=[]
    same_params=[]
    for var in params_data.variables:
        if params_data[var].dims != ('pft',):
            #print('Skipping var '+var)
            continue
        if var in ('mergetoclmpft','pftname','mxmat','pft','pftnum'):
            continue

        valdefault=params_default_data[var].values[pft_names_default.index(PFT_mappings[pftname])]
        valarctic=params_data[var].values[pft_names.index(pftname)]
        param_temptree=params_default_data[var].values[pft_names_default.index('broadleaf_evergreen_temperate_tree')]
        param_borshrub=params_default_data[var].values[pft_names_default.index('broadleaf_deciduous_boreal_shrub')]
        param_grass=params_default_data[var].values[pft_names_default.index('c3_non-arctic_grass')]
        param_arcgrass=params_default_data[var].values[pft_names_default.index('c3_arctic_grass')]

        if param_temptree != param_borshrub or param_grass != param_arcgrass:
            if hasattr(params_data[var],'units'):
                units=params_data[var].units
            elif hasattr(params_data[var],'unit'):
                units=params_data[var].unit
            else:
                units = '<Not provided>'
            if valdefault != valarctic and not (isnan(valarctic) and isnan(valdefault)):
                diff_params.append('Param {param:s} ({longname:s}, units = {units:s})'.format(param=var,longname=params_data[var].long_name,units=units))
                diff_params.append('   Default: {valdefault:1.3g}, Arctic: {valarctic:1.3g}'.format(valdefault=valdefault,valarctic=valarctic))
            else:
                same_params.append('Param {param:s} ({longname:s}), value = {val:1.3g}, units = {units:s}'.format(param=var,longname=params_data[var].long_name,val=valdefault,units=units))
    print('{nparams:d} parameters have different values:'.format(nparams=len(diff_params)//2))
    for line in diff_params:
        print(line)
    print('\n{nparams:d} parameters have the same values (but differ between other PFTs):'.format(nparams=len(same_params)//2))
    for line in same_params:
        print(line)

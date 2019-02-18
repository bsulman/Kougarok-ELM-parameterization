import xarray

basedir='/home/b0u/Kougarok_param_edits'
params=xarray.open_dataset(basedir+'/param_files/clm_params_newpfts_c180524_orig.nc',autoclose=True,)

pft_names=[name.strip() for name in params['pftname'].values.astype(str)]

landscape_ecotypes=['NAMC','DSLT','AS','WBT','TTWBT','TT']
ecotype_names={'NAMC':'Non-acidic mountain complex',
               'DSLT':'Dwarf shrub lichen tundra',
                'AS'  :'Alder shrubland',
                'WBT' :'Willow birch tundra',
                'TTWBT':'Tussock tundra/willow birch tundra',
                'TT'  :'Tussock tundra'}



def change_param(paramname,pftname,newval):
    print('** Changing parameter {paramname:s} for PFT {pftname:s} from {oldval:1.3g} to {newval:1.3g}'.format(
        oldval=params[paramname][pft_names.index(pftname)].values,
        newval=newval,paramname=paramname,pftname=pftname))
    params[paramname][pft_names.index(pftname)] = newval

def printnote(note):
    print('  *** Note: {note:s}'.format(note=note))

if __name__=='__main__':
    # Read in Verity Salmon's Kougarok measurements summary
    import pandas
    Koug_meas_biomass=pandas.read_excel(basedir+'/obs_data/NGEEArctic_Q3ELM_KougarokBiomass&NPP_20181112.xlsx',sheet_name='data')\
        .set_index(['Ecotype','ELMgroup'])
    Koug_meas_chem=pandas.read_excel(basedir+'/obs_data/NGEEArctic_Q3ELM_KougarokSLA&Chemistry_20181112.xlsx',sheet_name='data')\
        .rename(columns={'ELM_PFT':'ELMgroup'}).set_index(['Ecotype','ELMgroup'])

    surfdata=xarray.open_dataset(basedir+'/param_files/surfdata_51x63pt_kougarok-NGEE_TransA_simyr1850_c181115-sub12.nc')
    PFT_percents=pandas.DataFrame(data=surfdata.PCT_NAT_PFT.values.squeeze(),index=pft_names,columns=landscape_ecotypes)

    meas_leaf_C=(Koug_meas_biomass['LeafBiomass_gperm2']*Koug_meas_chem['LeafC_percent']/100)
    meas_root_C=(Koug_meas_biomass['FineRootBiomass_gperm2'][:,'mixed']*(Koug_meas_chem['FineRootC_percent']/100))[:,'mixed']

    obs_leafCN = Koug_meas_chem['LeafC_percent']/Koug_meas_chem['LeafN_percent']
    obs_stemCN = Koug_meas_chem['StemC_percent']/Koug_meas_chem['StemN_percent']
    obs_frootCN = 54.6/1.3 # Obs uses a single value for fine roots
    # Measured tissue C:N ratios are in weight units, and model expects gC/gN so they should be comparable.


    # For now: let's assume that relative amount of leaf biomass is proportional to relative amount of root biomass
    # But we may want to change to a different approach like PFT % coverage
    # Or does it make more sense to calculate these at the ecotype scale?
    # Also: This param is really defined to be NPP ratio, not biomass ratio. Should we parameterize using NPP measurements instead?
    leafCfrac=meas_leaf_C/meas_leaf_C.groupby('Ecotype').sum()

    def froot_leaf(ecotype,pft):
        froot_leaf=meas_root_C[ecotype]*leafCfrac[ecotype][pft]/meas_leaf_C[ecotype][pft]
        print('{ecotype:s}: {pft:s} leaf frac {leafCfrac:1.2f}, froot_leaf = {froot_leaf:1.2f}'.format(leafCfrac=leafCfrac[ecotype][pft],ecotype=ecotype,pft=pft,froot_leaf=froot_leaf))
        return froot_leaf

    # evergreen dwarf shrub
    froot_leaf_TT=froot_leaf('TT','dwarf shrub evergreen')
    froot_leaf_NAMC=froot_leaf('NAMC','dwarf shrub evergreen')
    change_param('froot_leaf','arctic_evergreen_shrub_dwarf',froot_leaf_NAMC )
    change_param('froot_long','arctic_evergreen_shrub_dwarf',1.5) # Value from Verity's data description Table 3

    # SLA in Verity's data is in cm2/g. Parameter in model is in m2/g. Divide obs by 100**2 to convert units
    change_param('slatop','arctic_evergreen_shrub_dwarf',Koug_meas_chem['LeafSLA_cm2perg'][:,'dwarf shrub evergreen'].mean()/100**2)
    change_param('leafcn','arctic_evergreen_shrub_dwarf',obs_leafCN[:,'dwarf shrub evergreen'].mean())
    change_param('frootcn','arctic_evergreen_shrub_dwarf',obs_frootCN)


    printnote('Model divides stems into "dead" (heartwood) and live components with different C:N ratios. How to compare with measurements?')

    # Dwarf deciduous shrub
    froot_leaf_DSLT=froot_leaf('DSLT','dwarf shrub deciduous')
    froot_leaf_WBT=froot_leaf('WBT','dwarf shrub deciduous')
    # Setting this to total root:leaf ratio of DSLT which is mostly shrubs.
    printnote('Setting dwarf deciduous shrub root_leaf to total root:leaf ratio for DSLT, which is mostly shrubs')
    change_param('froot_leaf','arctic_deciduous_shrub_dwarf',meas_root_C['DSLT']/meas_leaf_C['DSLT'].sum()  )

    change_param('slatop','arctic_deciduous_shrub_dwarf',Koug_meas_chem['LeafSLA_cm2perg'][:,'dwarf shrub deciduous'].mean()/100**2)
    change_param('leafcn','arctic_deciduous_shrub_dwarf',obs_leafCN[:,'dwarf shrub deciduous'].mean())
    change_param('frootcn','arctic_deciduous_shrub_dwarf',obs_frootCN)


    # Low deciduous shrub
    change_param('froot_leaf','arctic_deciduous_shrub_low',meas_root_C['DSLT']/meas_leaf_C['DSLT'].sum()  )

    change_param('slatop','arctic_deciduous_shrub_low',Koug_meas_chem['LeafSLA_cm2perg'][:,'low shrub deciduous'].mean()/100**2)
    change_param('leafcn','arctic_deciduous_shrub_low',obs_leafCN[:,'low shrub deciduous'].mean())
    change_param('frootcn','arctic_deciduous_shrub_low',obs_frootCN)


    # Tall non-alder shrub
    # Set this to total root_leaf ratio of all shrubs in WBT
    printnote('Using total root:leaf ratio of all shrubs in WBT for tall non-alder shrubs')
    leaf_shrubs_WBT = (meas_leaf_C['WBT'].sum()-meas_leaf_C['WBT']['graminoid'])
    froot_shrubs_WBT = meas_root_C['WBT']*(meas_leaf_C['WBT'].sum()-meas_leaf_C['WBT']['graminoid'])/meas_leaf_C['WBT'].sum()
    change_param('froot_leaf','arctic_deciduous_shrub_tall',froot_shrubs_WBT/leaf_shrubs_WBT)

    change_param('slatop','arctic_deciduous_shrub_tall',Koug_meas_chem['LeafSLA_cm2perg'].loc[:,['tall shrub deciduous willow','tall shrub deciduous birch']].mean()/100**2)
    printnote('Deciduous shrub measured C:N varies a lot, from 11 to 31. Mean is 22.')
    change_param('leafcn','arctic_deciduous_shrub_tall',obs_leafCN.loc[:,['tall shrub deciduous birch','tall shrub deciduous willow']].mean())
    change_param('frootcn','arctic_deciduous_shrub_tall',obs_frootCN)
    
    # Alder
    froot_leaf_AS=froot_leaf('AS','tall shrub deciduous alder')
    froot_leaf_TTWBT=froot_leaf('TTWBT','tall shrub deciduous alder')
    change_param('froot_leaf','arctic_deciduous_shrub_alder',froot_leaf_AS)

    change_param('slatop','arctic_deciduous_shrub_alder',Koug_meas_chem['LeafSLA_cm2perg'][:,'tall shrub deciduous alder'].mean()/100**2)
    change_param('leafcn','arctic_deciduous_shrub_alder',obs_leafCN[:,'tall shrub deciduous alder'].mean())
    change_param('frootcn','arctic_deciduous_shrub_alder',obs_frootCN)


    # Graminoid
    froot_leaf_TT=froot_leaf('TT','graminoid')
    froot_leaf_TTWBT=froot_leaf('TTWBT','graminoid')
    froot_leaf_WBT=froot_leaf('WBT','graminoid')
    # Use mean of TT and TTWBT, which have more graminoids and similar values
    # Use same values for wet and dry graminoids in model for now
    printnote('Using same parameter values for wet and dry graminoids')
    change_param('froot_leaf','arctic_dry_graminoid',0.5*(froot_leaf_TT+froot_leaf_TTWBT))
    change_param('froot_leaf','arctic_wet_graminoid',0.5*(froot_leaf_TT+froot_leaf_TTWBT))

    printnote('Graminoid SLA is much higher in WBT than other sites. What to do about that?')
    change_param('slatop','arctic_wet_graminoid',Koug_meas_chem['LeafSLA_cm2perg'][:,'graminoid'].mean()/100**2)
    change_param('slatop','arctic_dry_graminoid',Koug_meas_chem['LeafSLA_cm2perg'][:,'graminoid'].mean()/100**2)
    change_param('leafcn','arctic_dry_graminoid',obs_leafCN[:,'graminoid'].mean())
    change_param('leafcn','arctic_wet_graminoid',obs_leafCN[:,'graminoid'].mean())
    change_param('frootcn','arctic_dry_graminoid',obs_frootCN)
    change_param('frootcn','arctic_wet_graminoid',obs_frootCN)

    # Values from TT in Verity's data description Table 3
    change_param('froot_long','arctic_wet_graminoid',3.13)
    change_param('froot_long','arctic_dry_garminoid',3.13)

    # Forb
    # Probably not enough data to constrain roots (no site with high forb coverage). Make same as graminoids?
    # Leaving it alone for now.
    printnote('Not enough forb biomass at any site to estimate associated root biomass. Assume the ratio is the same as forbs?')
    printnote('No SLA measurements for forbs. Current forb value is {forbsla:1.2g}'.format(forbsla=params['slatop'].values[pft_names.index('arctic_forb')])) 

    
    change_param('leafcn','arctic_forb',obs_leafCN[:,'forb'].mean())
    change_param('frootcn','arctic_forb',obs_frootCN)

    print('Saving params file to clm_params_updated.nc')
    params.to_netcdf('clm_params_updated.nc')
    


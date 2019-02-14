import xarray

params=xarray.open_dataset('/lustre/or-hydra/cades-ccsi/proj-shared/project_acme/cesminput_ngee/lnd/clm2/paramdata/clm_params_c180524-sub12.nc',autoclose=True,)

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


if __name__=='__main__':
    # Read in Verity Salmon's Kougarok measurements summary
    import pandas
    Koug_meas_biomass=pandas.read_excel('/home/b0u/Kougarok_param_edits/NGEEArctic_Q3ELM_KougarokBiomass&NPP_20181112.xlsx',sheet_name='data')\
        .set_index(['Ecotype','ELMgroup'])
    Koug_meas_chem=pandas.read_excel('/home/b0u/Kougarok_param_edits/NGEEArctic_Q3ELM_KougarokSLA&Chemistry_20181112.xlsx',sheet_name='data')\
        .rename(columns={'ELM_PFT':'ELMgroup'}).set_index(['Ecotype','ELMgroup'])

    surfdata=xarray.open_dataset('/lustre/or-hydra/cades-ccsi/proj-shared/project_acme/cesminput_ngee/lnd/clm2/surfdata_map/surfdata_51x63pt_kougarok-NGEE_TransA_simyr1850_c181115-sub12.nc')
    PFT_percents=pandas.DataFrame(data=surfdata.PCT_NAT_PFT.values.squeeze(),index=pft_names,columns=landscape_ecotypes)

    meas_leaf_C=(Koug_meas_biomass['LeafBiomass_gperm2']*Koug_meas_chem['LeafC_percent']/100)
    meas_root_C=(Koug_meas_biomass['FineRootBiomass_gperm2'][:,'mixed']*(Koug_meas_chem['FineRootC_percent']/100))[:,'mixed']

    # For now: let's assume that relative amount of leaf biomass is proportional to relative amount of root biomass
    # But we may want to change to a different approach like PFT % coverage
    leafCfrac=meas_leaf_C/meas_leaf_C.groupby('Ecotype').sum()

    def froot_leaf(ecotype,pft):
        froot_leaf=meas_root_C[ecotype]*leafCfrac[ecotype][pft]/meas_leaf_C[ecotype][pft]
        print('{ecotype:s}: {pft:s} leaf frac {leafCfrac:1.2f}, froot_leaf = {froot_leaf:1.2f}'.format(leafCfrac=leafCfrac[ecotype][pft],ecotype=ecotype,pft=pft,froot_leaf=froot_leaf))
        return froot_leaf

    # evergreen dwarf shrub
    froot_leaf_TT=froot_leaf('TT','dwarf shrub evergreen')
    froot_leaf_NAMC=froot_leaf('NAMC','dwarf shrub evergreen')
    change_param('froot_leaf','arctic_evergreen_shrub_dwarf',froot_leaf_NAMC )

    # Alder
    froot_leaf_AS=froot_leaf('AS','tall shrub deciduous alder')
    froot_leaf_TTWBT=froot_leaf('TTWBT','tall shrub deciduous alder')
    change_param('froot_leaf','arctic_deciduous_shrub_alder',froot_leaf_AS)






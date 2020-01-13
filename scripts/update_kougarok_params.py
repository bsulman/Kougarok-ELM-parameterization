import xarray

basedir='..'
params=xarray.open_dataset(basedir+'/param_files/clm_params_newpfts_c180524_orig.nc')
params=xarray.open_dataset(basedir+'/param_files/clm_params_defaultpfts_c180524_orig.nc') 

pft_names=[name.strip() for name in params['pftname'].values.astype(str)]

landscape_ecotypes=['NAMC','DSLT','AS','WBT','TTWBT','TT']
ecotype_names={'NAMC':'Non-acidic mountain complex',
               'DSLT':'Dwarf shrub lichen tundra',
                'AS'  :'Alder shrubland',
                'WBT' :'Willow birch tundra',
                'TTWBT':'Tussock tundra/willow birch tundra',
                'TT'  :'Tussock tundra'}

def pft_params(paramdata,paramnames):
    pftnames=[name.strip() for name in paramdata['pftname'].values.astype(str)]
    if isinstance(paramnames,str):
        paramnames=[paramnames]
    pdict={'pftname':pftnames}
    for name in paramnames:
       pdict[name]=paramdata[name]
    return pandas.DataFrame(pdict).set_index('pftname')

def change_param(paramname,pftname,newval):
    if 'pft' not in params[paramname].coords:
        raise RuntimeError('Parameter {param:s} does not have a PFT dimension. Try using change_universal_param instead'.format(param=paramname))
    print('** Changing parameter {paramname:s} for PFT {pftname:s} from {oldval:1.3g} to {newval:1.3g}'.format(
        oldval=params[paramname][pft_names.index(pftname)].values,
        newval=newval,paramname=paramname,pftname=pftname))
    params[paramname][pft_names.index(pftname)] = newval

def change_universal_param(paramname,newval):
    if 'pft' in params[paramname].coords:
        raise RuntimeError('Parameter {param:s} has a PFT dimension. Are you sure you want to change all of them?'.format(param=paramname))
    print('** Changing parameter {paramname:s} for ALL PFTS from {oldval:1.3g} to {newval:1.3g}'.format(
        oldval=params[paramname].values[0],
        newval=newval,paramname=paramname))
    params[paramname] = newval

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
    meas_rhizome_C=(Koug_meas_biomass['RhizomeBiomass_gperm2']*Koug_meas_chem['RhizomeC_percent']/100)
    meas_rhizome_NPP=(Koug_meas_biomass['RhizomeNPP_gperm2peryr']*Koug_meas_chem['RhizomeC_percent']/100)
    meas_leaf_NPP   =(Koug_meas_biomass['LeafNPP_gperm2peryr']*Koug_meas_chem['LeafC_percent']/100)
    meas_stem_NPP   =(Koug_meas_biomass['StemNPP_gperm2peryr']*Koug_meas_chem['StemC_percent']/100)
    meas_root_NPP   =(Koug_meas_biomass['FineRootNPP_gperm2peryr']*Koug_meas_chem['FineRootC_percent']/100)[:,'mixed']

    obs_leafCN = Koug_meas_chem['LeafC_percent']/Koug_meas_chem['LeafN_percent']
    obs_stemCN = Koug_meas_chem['StemC_percent']/Koug_meas_chem['StemN_percent']
    obs_frootCN = 58.0 # Mean of measured root C:N across Kougarok plots
    obs_rhizomeCN = Koug_meas_chem['RhizomeC_percent']/Koug_meas_chem['RhizomeN_percent']
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


    meas_tot_NPP    = meas_rhizome_NPP+meas_leaf_NPP+meas_stem_NPP+meas_root_NPP*leafCfrac

    # fcur_deciduous=0.0
    # fcur_evergreen=1.0
    printnote('Setting fcur to one minus the rhizome fraction of total NPP')
    fcur_pfts = 1.0 - meas_rhizome_NPP/meas_tot_NPP
    fcur_pfts[:]=0.5

    # evergreen dwarf shrub
    pft='arctic_evergreen_shrub_dwarf'
    froot_leaf_TT=froot_leaf('TT','dwarf shrub evergreen')
    froot_leaf_NAMC=froot_leaf('NAMC','dwarf shrub evergreen')
    # printnote('Setting froot_leaf for evergreen species to include rhizomes, and weighting root turnover by rhizome biomass and turnover rates')
    printnote('Treating rhizomes as ELM storage C and N pool, not as roots')
    printnote('Assuming a rhizome (storage) turnover of 5 years based on dwarf evergreen shrub value from Table 5')
    printnote('In the future, this should be a PFT-specific parameter')
    change_universal_param('fstor2tran',1.0/4.0)

    leaflong=(Koug_meas_biomass['LeafBiomass_gperm2']/Koug_meas_biomass['LeafNPP_gperm2peryr'])[:,'dwarf shrub evergreen'].mean()
    leaflong=3.5
    change_param('leaf_long',pft,leaflong )


    # Dwarf evergreen shrub rhizome turnover time of 5 years from Table 5 in Verity's data description
    # rhizome_leaf_NAMC=meas_rhizome_C['NAMC']['dwarf shrub evergreen']/meas_leaf_C['NAMC']['dwarf shrub evergreen']
    # dwarf_e_shrub_frootlong=(1.5*froot_leaf_NAMC + 20.0*rhizome_leaf_NAMC)/(froot_leaf_NAMC+rhizome_leaf_NAMC)
    dwarf_e_shrub_frootlong=1.56 # Table 3 in Verity's metadata
    dwarf_e_shrub_frootlong=2.0
    change_param('froot_long',pft,dwarf_e_shrub_frootlong) # Longevity estimated from Verity's data description Table 3
    # change_param('froot_leaf',pft,(froot_leaf_NAMC)*leaflong/dwarf_e_shrub_frootlong )
    change_param('froot_leaf',pft,3.0)
    # fcur will have to be calibrated so storage pool is consistent with measured rhizome biomass
    change_param('fcur',pft,fcur_pfts['NAMC','dwarf shrub evergreen'])

    # SLA in Verity's data is in cm2/g. Parameter in model is in m2/g. Divide obs by 100**2 to convert units
    change_param('slatop',pft,Koug_meas_chem['LeafSLA_cm2perg'][:,'dwarf shrub evergreen'].mean()/100**2)
    change_param('leafcn',pft,obs_leafCN[:,'dwarf shrub evergreen'].mean())
    # change_param('frootcn',pft,(obs_frootCN*froot_leaf_NAMC + obs_rhizomeCN[:,'dwarf shrub evergreen'].mean()*rhizome_leaf_NAMC)/(froot_leaf_NAMC+rhizome_leaf_NAMC))
    change_param('frootcn',pft,obs_frootCN)
    change_param('croot_stem',pft,0.1)
    change_param('stem_leaf',pft,0.1)

    printnote('Model divides stems into "dead" (heartwood) and live components with different C:N ratios. How to compare with measurements?')

    # Dwarf deciduous shrub
    pft='arctic_deciduous_shrub_dwarf'
    froot_leaf_DSLT=froot_leaf('DSLT','dwarf shrub deciduous')
    froot_leaf_WBT=froot_leaf('WBT','dwarf shrub deciduous')
    # Setting this to total root:leaf ratio of DSLT which is mostly shrubs.
    printnote('Setting dwarf deciduous shrub root_leaf to total root:leaf ratio for DSLT, which is mostly shrubs')
    # This should be adjusted to reflect longer lifetime of fine roots
    change_param('froot_leaf',pft,meas_root_C['DSLT']/meas_leaf_C['DSLT'].sum()/dwarf_e_shrub_frootlong  )
    change_param('froot_long',pft,dwarf_e_shrub_frootlong) # Not sure if this has any effect for deciduous species

    change_param('slatop',pft,Koug_meas_chem['LeafSLA_cm2perg'][:,'dwarf shrub deciduous'].mean()/100**2)
    change_param('leafcn',pft,obs_leafCN[:,'dwarf shrub deciduous'].mean())
    change_param('frootcn',pft,obs_frootCN)
    change_param('fcur',pft,fcur_pfts['DSLT','dwarf shrub deciduous'])

    # Tall non-alder shrub
    pft='arctic_deciduous_shrub_tall'
    # Set this to total root_leaf ratio of all shrubs in WBT
    printnote('Using total root:leaf ratio of all shrubs in WBT for tall non-alder shrubs')
    leaf_shrubs_WBT = (meas_leaf_C['WBT'].sum()-meas_leaf_C['WBT']['graminoid'])
    froot_shrubs_WBT = meas_root_C['WBT']*(meas_leaf_C['WBT'].sum()-meas_leaf_C['WBT']['graminoid'])/meas_leaf_C['WBT'].sum()
    change_param('froot_leaf',pft,froot_shrubs_WBT/leaf_shrubs_WBT)

    change_param('slatop',pft,Koug_meas_chem['LeafSLA_cm2perg'].loc[:,['tall shrub deciduous willow','tall shrub deciduous birch']].mean()/100**2)
    printnote('Deciduous shrub measured C:N varies a lot, from 11 to 31. Mean is 22.')
    change_param('leafcn',pft,obs_leafCN.loc[:,['tall shrub deciduous birch','tall shrub deciduous willow']].mean())
    change_param('frootcn',pft,obs_frootCN)
    change_param('fcur',pft,fcur_pfts[:,'tall shrub deciduous willow'].mean())
    
    
    # Low deciduous shrub
    pft='arctic_deciduous_shrub_low'
    #change_param('froot_leaf',pft,meas_root_C['DSLT']/meas_leaf_C['DSLT'].sum()  )
    printnote('Making froot_leaf for low shrubs intermediate between dwarf and tall values')
    change_param('froot_leaf',pft, 0.5*(params['froot_leaf'][pft_names.index('arctic_deciduous_shrub_dwarf')]+params['froot_leaf'][pft_names.index('arctic_deciduous_shrub_tall')]).values)
    change_param('froot_long',pft,dwarf_e_shrub_frootlong)

    change_param('slatop',pft,Koug_meas_chem['LeafSLA_cm2perg'][:,'low shrub deciduous'].mean()/100**2)
    change_param('leafcn',pft,obs_leafCN[:,'low shrub deciduous'].mean())
    change_param('frootcn',pft,obs_frootCN)
    change_param('fcur',pft,fcur_pfts[:,'low shrub deciduous'].mean() )

    
    # Alder
    pft='arctic_deciduous_shrub_alder'
    froot_leaf_AS=froot_leaf('AS','tall shrub deciduous alder')
    froot_leaf_TTWBT=froot_leaf('TTWBT','tall shrub deciduous alder')
    change_param('froot_leaf',pft,froot_leaf_AS)

    change_param('slatop',pft,Koug_meas_chem['LeafSLA_cm2perg'][:,'tall shrub deciduous alder'].mean()/100**2)
    change_param('leafcn',pft,obs_leafCN[:,'tall shrub deciduous alder'].mean())
    change_param('frootcn',pft,obs_frootCN)
    change_param('fcur',pft,fcur_pfts[:,'tall shrub deciduous alder'].mean() )
    change_param('stem_leaf',pft,0.20)
    change_param('croot_stem',pft,0.5)


    # Graminoid
    pftwet='arctic_wet_graminoid'
    pftdry='arctic_dry_graminoid'
    froot_leaf_TT=froot_leaf('TT','graminoid')
    froot_leaf_TTWBT=froot_leaf('TTWBT','graminoid')
    froot_leaf_WBT=froot_leaf('WBT','graminoid')
    # rhizome_leaf_TT=meas_rhizome_C['TT']['graminoid']/meas_leaf_C['TT']['graminoid']
    # rhizome_leaf_TTWBT=meas_rhizome_C['TTWBT']['graminoid']/meas_leaf_C['TTWBT']['graminoid']

    # Use mean of TT and TTWBT, which have more graminoids and similar values
    # Use same values for wet and dry graminoids in model for now
    printnote('Using same parameter values for wet and dry graminoids')
    # printnote('Calculating rhizome lifetime from NPP and biomass (assuming steady state)')
    # rhizome_lifetime = (meas_rhizome_C/meas_rhizome_NPP)['TT','graminoid']
    # froot_leaf_gram=0.5*(froot_leaf_TT+froot_leaf_TTWBT)/3.13 + 0.5*(rhizome_leaf_TT+rhizome_leaf_TTWBT)/rhizome_lifetime
    # froot_leaf_gram=froot_leaf_gram*1.5
    froot_leaf_gram = 0.5*froot_leaf_TT + 0.5*froot_leaf_TTWBT

    # Numbers for root and leaf longevity
    printnote('Assigning graminoid leaf longevity of 2 years based on 50% leaf replacement estimate from Shaver and Laundre 2003 GCB paper')
    # Root longevity values from TT in Verity's data description Table 3
    rootlong=3.13
    leaflong=2.0
    change_param('froot_leaf',pftdry,froot_leaf_gram*leaflong/rootlong*1.2)
    change_param('froot_leaf',pftwet,froot_leaf_gram*leaflong/rootlong*1.2)

    printnote('Graminoid SLA is much higher in WBT than other sites. What to do about that?')
    change_param('slatop',pftwet,Koug_meas_chem['LeafSLA_cm2perg'][:,'graminoid'].mean()/100**2)
    change_param('slatop',pftdry,Koug_meas_chem['LeafSLA_cm2perg'][:,'graminoid'].mean()/100**2)
    change_param('leafcn',pftdry,obs_leafCN[:,'graminoid'].mean())
    change_param('leafcn',pftwet,obs_leafCN[:,'graminoid'].mean())
    printnote('Setting graminoid root C:N based on updated Kougarok chemistry')
    change_param('frootcn',pftdry,75.0)
    change_param('frootcn',pftwet,75.0)

    #frootfrac_graminoid=(froot_leaf_TT+froot_leaf_TTWBT)/(rhizome_leaf_TT+rhizome_leaf_TTWBT+froot_leaf_TT+froot_leaf_TTWBT)
    #change_param('froot_long',pftwet,3.13*frootfrac_graminoid + rhizome_lifetime*(1-frootfrac_graminoid))
    #change_param('froot_long',pftdry,3.13*frootfrac_graminoid + rhizome_lifetime*(1-frootfrac_graminoid))
    change_param('froot_long',pftwet,rootlong)
    change_param('froot_long',pftdry,rootlong)

    change_param('fcur',pftwet,fcur_pfts[:,'graminoid'].mean() )
    change_param('fcur',pftdry,fcur_pfts[:,'graminoid'].mean() )

    printnote('According to Verity, grasses act more like evergreen plants. Do not drop leaves/roots every year')
    change_param('season_decid',pftwet,0)
    change_param('season_decid',pftdry,0)
    change_param('evergreen',pftwet,1)
    change_param('evergreen',pftdry,1)

    change_param('leaf_long',pftwet,leaflong)
    change_param('leaf_long',pftdry,leaflong)

    printnote('VCmax for graminoids seems too high with updated parameters. Reducing it by about 50%')
    change_param('flnr',pftwet,0.09)
    change_param('flnr',pftdry,0.09)

    # Forb
    # Probably not enough data to constrain roots (no site with high forb coverage). Make same as graminoids?
    # Leaving it alone for now.
    printnote('Not enough forb biomass at any site to estimate associated root biomass. Assume the ratio is the same as grasses?')
    printnote('No SLA measurements for forbs. Current forb value is {forbsla:1.2g}'.format(forbsla=params['slatop'].values[pft_names.index('arctic_forb')])) 
    printnote('Should forbs be evergreen or deciduous?')

    pft='arctic_forb' 
    change_param('leafcn',pft,obs_leafCN[:,'forb'].mean())
    change_param('frootcn',pft,obs_frootCN)
 
    change_param('fcur',pft,fcur_pfts[:,'forb'].mean())
    change_param('froot_leaf',pft,2.0)
    change_param('flnr',pft,0.2)

    # Lichen
    change_param('froot_leaf','arctic_lichen',0.2)
    change_param('leaf_long','arctic_lichen',10.0)

    # Set up params for new dormant maintenance respiration
    dormant_mr_temp=273.15+2.5
    dormant_mr_factor=5e-2
    printnote('Setting dormancy temperature to {0:1.1f} C'.format(dormant_mr_temp-273.15))
    printnote('Setting dormancy maintenance resp factor to {0:1.1g}'.format(dormant_mr_factor))
    params['dormant_mr_temp']=xarray.DataArray(name='dormant_mr_temp',dims='allpfts',data=[dormant_mr_temp],attrs={'units':'degrees K','long_name':'Maximum temperature for dormant maintenance respiration'})
    params['dormant_mr_factor']=xarray.DataArray(name='dormant_mr_factor',dims='allpfts',data=[dormant_mr_factor],attrs={'units':'unitless','long_name':'Dormant maintenance respiration multiplication factor'})

    # Set up N fixation params
    from numpy import zeros
    params['Nfix_NPP_c1']=xarray.DataArray(name='Nfix_NPP_c1',dims='pft',data=zeros(len(pft_names))+1.8,attrs={'units':'gN/m2/yr','long_name':'Pre-exponential factor in NPP-Nfix equation'})
    params['Nfix_NPP_c2']=xarray.DataArray(name='Nfix_NPP_c2',dims='pft',data=zeros(len(pft_names))+0.003,attrs={'units':'gN/m2/yr','long_name':'Exponential factor in NPP-Nfix equation'})
    change_param('Nfix_NPP_c1','arctic_deciduous_shrub_alder',10.0)
    change_param('Nfix_NPP_c1','arctic_forb',5.0)
    change_param('Nfix_NPP_c2','arctic_deciduous_shrub_alder',0.05)

    print('Saving params file to clm_params_updated.nc')
    params.to_netcdf(basedir+'/clm_params_updated.nc',format='NETCDF4_CLASSIC')
    



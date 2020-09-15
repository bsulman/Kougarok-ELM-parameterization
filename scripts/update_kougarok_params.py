import xarray
import numpy
from numpy import zeros

basedir='..'
params_newfengming=xarray.open_dataset(basedir+'/param_files/clm_params_newpfts_c180524_orig.nc',decode_times=False)
params_orig=xarray.open_dataset(basedir+'/param_files/clm_params_defaultpfts_c180524_orig.nc',decode_times=False) 

reset_to_defaults=True

pft_names=[name.strip() for name in params_newfengming['pftname'].values.astype(str)]

landscape_ecotypes=['NAMC','DSLT','AS','WBT','TTWBT','TT']
ecotype_names={'NAMC':'Non-acidic mountain complex',
               'DSLT':'Dwarf shrub lichen tundra',
                'AS'  :'Alder shrubland',
                'WBT' :'Willow birch tundra',
                'TTWBT':'Tussock tundra/willow birch tundra',
                'TT'  :'Tussock tundra'}




new_old_mapping={
    'not_vegetated':'not_vegetated',
    'arctic_lichen':'c3_arctic_grass',
    'arctic_bryophyte':'c3_arctic_grass',
    'arctic_evergreen_shrub_dwarf':'broadleaf_deciduous_boreal_shrub',
    'arctic_evergreen_shrub_tall':'broadleaf_deciduous_boreal_shrub',
    'arctic_deciduous_shrub_dwarf':'broadleaf_deciduous_boreal_shrub',
    'arctic_deciduous_shrub_low':'broadleaf_deciduous_boreal_shrub',
    'arctic_deciduous_shrub_tall':'broadleaf_deciduous_boreal_shrub',
    'arctic_deciduous_shrub_alder':'broadleaf_deciduous_boreal_shrub',
    'arctic_forb':'c3_arctic_grass',
    'arctic_dry_graminoid':'c3_arctic_grass',
    'arctic_wet_graminoid':'c3_arctic_grass',
}
names_new=[s.strip() for s in params_newfengming['pftname'].values.astype(str)]
names_old=[s.strip() for s in params_orig['pftname'].values.astype(str)]
mapping_array=numpy.array([names_old.index(new_old_mapping[n]) for n in names_new])

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

    Koug_meas_biomass_new=pandas.read_excel(basedir+'/obs_data/NGEE Arctic Veg data compiled 12June2019/NGEEArctic_Q3ELM_KougarokBiomass&NPP_20190612.xlsx',sheet_name='data')\
        .set_index(['Ecotype','PlotID','ELM_PFT'])
    Koug_meas_chem_new=pandas.read_excel(basedir+'/obs_data/NGEE Arctic Veg data compiled 12June2019/Kougarok_Q3ELM_KougarokSLA&Chemistry_20190612.xlsx',sheet_name='data')\
        .set_index(['Ecotype','PlotID','ELM_PFT'])


    surfdata=xarray.open_dataset(basedir+'/param_files/surfdata_51x63pt_kougarok-NGEE_TransA_simyr1850_c181115-sub12.nc')
    PFT_percents=pandas.DataFrame(data=surfdata.PCT_NAT_PFT.values.squeeze(),index=pft_names,columns=landscape_ecotypes)

    params=params_newfengming.copy()
    if reset_to_defaults:
        # Reset all values to equivalent original ones
        for var in params.variables:
            if var in ['mergetoclmpft','pft','pftname','pftnum']:
                continue
            if params[var].dims == ('pft',):
                params[var]=params_orig[var][mapping_array]
            elif len(params[var].dims)>0 and params[var].dims[0] == 'pft' and len(params[var].dims) == 2: # Should only be one thing that is (pft,levgrnd)
                params[var]=params_orig[var][mapping_array,:]
            elif 'pft' in params[var].dims:
                raise ValueError('No rule for copying dimension of this type: %s'%str(params[var].dims))
            else:
                # Just copy it over since it is not pft-specific
                params[var]=params_orig[var]

        for pftnum in range(len(params['pft'])):
            pftname=names_new[pftnum]
            if 'evergreen_shrub' in pftname or 'lichen' in pftname or 'bryophyte' in pftname:
                change_param('season_decid',pftname,0)
                change_param('evergreen',pftname,1)
                change_param('fcur',pftname,1.0)
                change_param('flnr',pftname,0.0755)
                change_param('pftpar30',pftname,0.0)
        
        change_param('froot_leaf','arctic_lichen',0.1)
        change_param('froot_leaf','arctic_bryophyte',0.1)
        change_param('flnr','arctic_lichen',0.0435)
        change_param('flnr','arctic_bryophyte',0.0485)


    params['rhizome_long']=xarray.DataArray(name='rhizome_long',dims='pft',data=zeros(len(pft_names))+2.0,attrs={'units':'yr','long_name':'Nonwoody rhizome longevity'})

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
    # SLA in Verity's data is in cm2/g. Parameter in model is in m2/gC. Divide obs by 100**2 and %C to convert units
    meas_SLA = Koug_meas_chem['LeafSLA_cm2perg']/(Koug_meas_chem['LeafC_percent']/100)/100**2

    fcur_deciduous=0.0
    fcur_evergreen=0.75

    # evergreen dwarf shrub
    pft='arctic_evergreen_shrub_dwarf'
    froot_leaf_TT=froot_leaf('TT','dwarf shrub evergreen')
    froot_leaf_NAMC=froot_leaf('NAMC','dwarf shrub evergreen')
    # printnote('Setting froot_leaf for evergreen species to include rhizomes, and weighting root turnover by rhizome biomass and turnover rates')
    # printnote('Treating rhizomes as ELM storage C and N pool, not as roots')
    # printnote('Assuming a rhizome (storage) turnover of 5 years based on dwarf evergreen shrub value from Table 5')
    # printnote('In the future, this should be a PFT-specific parameter')
    change_universal_param('fstor2tran',1.0/4.0)

    leaflong=(Koug_meas_biomass['LeafBiomass_gperm2']/Koug_meas_biomass['LeafNPP_gperm2peryr'])[:,'dwarf shrub evergreen'].mean()
    leaflong=3.5
    change_param('leaf_long',pft,leaflong )
    change_param('s_vc',pft,20.72) # Consistent with needleleaf trees
    change_param('leafcn_obs',pft,22.0)
    change_param('leafcn_obs','arctic_evergreen_shrub_tall',22.0)

    change_param('rholnir',pft,0.35)
    change_param('rholnir','arctic_evergreen_shrub_tall',0.35)
    change_param('taulnir',pft,0.10)
    change_param('taulnir','arctic_evergreen_shrub_tall',0.10)
    change_param('rholvis',pft,0.07)
    change_param('rholvis','arctic_evergreen_shrub_tall',0.07)

    # Dwarf evergreen shrub rhizome turnover time of 5 years from Table 5 in Verity's data description
    # rhizome_leaf_NAMC=meas_rhizome_C['NAMC']['dwarf shrub evergreen']/meas_leaf_C['NAMC']['dwarf shrub evergreen']
    # dwarf_e_shrub_frootlong=(1.5*froot_leaf_NAMC + 20.0*rhizome_leaf_NAMC)/(froot_leaf_NAMC+rhizome_leaf_NAMC)
    dwarf_e_shrub_frootlong=1.56 # Table 3 in Verity's metadata
    dwarf_e_shrub_frootlong=2.0
    change_param('froot_long',pft,dwarf_e_shrub_frootlong) # Longevity estimated from Verity's data description Table 3
    # change_param('froot_leaf',pft,(froot_leaf_NAMC)*leaflong/dwarf_e_shrub_frootlong )
    change_param('froot_leaf',pft,2.5)
    change_param('fcur',pft,fcur_evergreen)

    change_param('slatop',pft,meas_SLA[:,'dwarf shrub evergreen'].mean())
    change_param('leafcn',pft,obs_leafCN[:,'dwarf shrub evergreen'].mean())
    # change_param('frootcn',pft,(obs_frootCN*froot_leaf_NAMC + obs_rhizomeCN[:,'dwarf shrub evergreen'].mean()*rhizome_leaf_NAMC)/(froot_leaf_NAMC+rhizome_leaf_NAMC))
    change_param('frootcn',pft,obs_frootCN)
    change_param('croot_stem',pft,1.0)
    change_param('stem_leaf',pft,0.15)

    change_param('lflitcn',pft,56.0)

    printnote('Model divides stems into "dead" (heartwood) and live components with different C:N ratios. How to compare with measurements?')

    # Dwarf deciduous shrub
    pft='arctic_deciduous_shrub_dwarf'
    froot_leaf_DSLT=froot_leaf('DSLT','dwarf shrub deciduous')
    froot_leaf_WBT=froot_leaf('WBT','dwarf shrub deciduous')
    # Setting this to total root:leaf ratio of DSLT which is mostly shrubs.
    printnote('Setting dwarf deciduous shrub root_leaf to total root:leaf ratio for DSLT, which is mostly shrubs')
    # This should be adjusted to reflect longer lifetime of fine roots
    #change_param('froot_leaf',pft,meas_root_C['DSLT']/meas_leaf_C['DSLT'].sum()/dwarf_e_shrub_frootlong  )
    change_param('froot_leaf',pft,1.5)
    change_param('froot_long',pft,dwarf_e_shrub_frootlong) # Not sure if this has any effect for deciduous species
    change_param('croot_stem',pft,0.6)
    change_param('stem_leaf',pft,0.22)

    change_param('slatop',pft,meas_SLA[:,'dwarf shrub deciduous'].mean())
    change_param('leafcn',pft,obs_leafCN[:,'dwarf shrub deciduous'].mean())
    change_param('frootcn',pft,obs_frootCN)
    change_param('fcur',pft,0.0)

    # Tall non-alder shrub
    pft='arctic_deciduous_shrub_tall'
    # Set this to total root_leaf ratio of all shrubs in WBT
    printnote('Using total root:leaf ratio of all shrubs in WBT for tall non-alder shrubs')
    leaf_shrubs_WBT = (meas_leaf_C['WBT'].sum()-meas_leaf_C['WBT']['graminoid'])
    froot_shrubs_WBT = meas_root_C['WBT']*(meas_leaf_C['WBT'].sum()-meas_leaf_C['WBT']['graminoid'])/meas_leaf_C['WBT'].sum()
    change_param('froot_leaf',pft,froot_shrubs_WBT/leaf_shrubs_WBT)

    change_param('slatop',pft,meas_SLA.loc[:,['tall shrub deciduous willow','tall shrub deciduous birch']].mean())
    printnote('Deciduous shrub measured C:N varies a lot, from 11 to 31. Mean is 22.')
    change_param('leafcn',pft,obs_leafCN.loc[:,['tall shrub deciduous birch','tall shrub deciduous willow']].mean())
    change_param('frootcn',pft,obs_frootCN)
    change_param('fcur',pft,fcur_deciduous)
    change_param('croot_stem',pft,0.5)
    
    
    # Low deciduous shrub
    pft='arctic_deciduous_shrub_low'
    #change_param('froot_leaf',pft,meas_root_C['DSLT']/meas_leaf_C['DSLT'].sum()  )
    printnote('Making froot_leaf for low shrubs intermediate between dwarf and tall values')
    #change_param('froot_leaf',pft, 0.5*(params['froot_leaf'][pft_names.index('arctic_deciduous_shrub_dwarf')]+params['froot_leaf'][pft_names.index('arctic_deciduous_shrub_tall')]).values)
    change_param('froot_leaf',pft,1.41)
    change_param('froot_long',pft,dwarf_e_shrub_frootlong)

    change_param('slatop',pft,meas_SLA[:,'low shrub deciduous'].mean())
    change_param('leafcn',pft,obs_leafCN[:,'low shrub deciduous'].mean())
    change_param('frootcn',pft,obs_frootCN)
    change_param('fcur',pft,fcur_deciduous )
    change_param('croot_stem',pft,0.8)
    change_param('stem_leaf',pft,0.22)
    
    # Alder
    pft='arctic_deciduous_shrub_alder'
    froot_leaf_AS=froot_leaf('AS','tall shrub deciduous alder')
    froot_leaf_TTWBT=froot_leaf('TTWBT','tall shrub deciduous alder')
    change_param('froot_leaf',pft,0.25)

    change_param('slatop',pft,meas_SLA[:,'tall shrub deciduous alder'].mean())
    change_param('leafcn',pft,obs_leafCN[:,'tall shrub deciduous alder'].mean())
    change_param('frootcn',pft,obs_frootCN)
    change_param('fcur',pft,fcur_deciduous )
    change_param('stem_leaf',pft,0.5)
    change_param('croot_stem',pft,0.6)
    # change_param('flnr',pft,0.1365*2) # Set to double of other shrubs
    # Probably not the solution... Alistair's measurements show alder vcmax of ~50-100, which is consistent
    # with model calculation using default flnr value. Doubled flnr gives vcmax of 200, which seems too high

    # Graminoid
    froot_leaf_TT=froot_leaf('TT','graminoid')
    froot_leaf_TTWBT=froot_leaf('TTWBT','graminoid')
    froot_leaf_WBT=froot_leaf('WBT','graminoid')
    # rhizome_leaf_TT=meas_rhizome_C['TT']['graminoid']/meas_leaf_C['TT']['graminoid']
    # rhizome_leaf_TTWBT=meas_rhizome_C['TTWBT']['graminoid']/meas_leaf_C['TTWBT']['graminoid']

    # Use mean of TT and TTWBT, which have more graminoids and similar values
    # Use same values for wet and dry graminoids in model for now
    printnote('Using same parameter values for wet and dry graminoids')
    printnote('Calculating rhizome lifetime from NPP and biomass (assuming steady state)')
    rhizome_lifetime = (meas_rhizome_C/meas_rhizome_NPP)['TT','graminoid']
    # froot_leaf_gram=0.5*(froot_leaf_TT+froot_leaf_TTWBT)/3.13 + 0.5*(rhizome_leaf_TT+rhizome_leaf_TTWBT)/rhizome_lifetime
    # froot_leaf_gram=froot_leaf_gram*1.5
    froot_leaf_gram = 0.5*froot_leaf_TT + 0.5*froot_leaf_TTWBT

    # Numbers for root and leaf longevity
    printnote('Assigning graminoid leaf longevity of 2 years based on 50% leaf replacement estimate from Shaver and Laundre 2003 GCB paper')
    # Root longevity values from TT in Verity's data description Table 3
    rootlong=3.13
    leaflong=2.0
    # change_param('froot_leaf',pftdry,froot_leaf_gram*leaflong/rootlong*1.2)
    # change_param('froot_leaf',pftwet,froot_leaf_gram*leaflong/rootlong*1.2)
    params['lwtop_pfts']=xarray.DataArray(name='lwtop_pfts',dims='pft',data=zeros(len(pft_names))+params['lwtop_ann'].values,attrs={'units':'unitless','long_name':'Live wood turnover proportion (PFT-specific)'})
    for pft in ['arctic_wet_graminoid','arctic_dry_graminoid']:
        change_param('froot_leaf',pft,4.0)

        # printnote('Graminoid SLA is much higher in WBT than other sites. What to do about that?')
        change_param('slatop',pft,meas_SLA[:,'graminoid'].mean())
        change_param('leafcn',pft,obs_leafCN[:,'graminoid'].mean())
        printnote('Setting graminoid root C:N based on updated Kougarok chemistry')
        change_param('frootcn',pft,75.0)

        #frootfrac_graminoid=(froot_leaf_TT+froot_leaf_TTWBT)/(rhizome_leaf_TT+rhizome_leaf_TTWBT+froot_leaf_TT+froot_leaf_TTWBT)
        #change_param('froot_long',pftwet,3.13*frootfrac_graminoid + rhizome_lifetime*(1-frootfrac_graminoid))
        #change_param('froot_long',pftdry,3.13*frootfrac_graminoid + rhizome_lifetime*(1-frootfrac_graminoid))
        change_param('froot_long',pft,rootlong)
        change_param('rhizome_long',pft,rhizome_lifetime)

        change_param('fcur',pft,fcur_evergreen )

        printnote('According to Verity, grasses act more like evergreen plants. Do not drop leaves/roots every year')
        change_param('season_decid',pft,0)
        change_param('evergreen',pft,1)

        change_param('leaf_long',pft,leaflong)

        printnote('VCmax for graminoids seems too high with updated parameters. Reducing it by about 50%')
        change_param('flnr',pft,0.09)

        ## printnote('Making graminoids into woody plants to allow rhizomes. Setting with very low stem and very high live wood frac')
        ## Update: this doesn't work. Making them woody plants without deadwood stems messes up height calculations etc
        ## change_param('woody',pft,1)
        ## change_param('croot_stem',pft,10.0)  # Actual allocation to Croots is leaf_alloc*stem_leaf*croot_stem
        ## change_param('flivewd',pft,0.99) # Live wood fraction
        ## change_param('lwtop_pfts',pft,0.0) # No turnover to dead wood for graminoids    
        ## change_param('deadwdcn',pft,500)   # If this is zero the model crashes on initialization

        # Updated strategy: Changed code so for non-woody plants stem_leaf is interpreted as rhizome allocation and sent to livecrootc 
        change_param('stem_leaf',pft,0.2)
        change_param('livewdcn',pft,50)   # If this is zero the model crashes on initialization

    # Forb
    # Probably not enough data to constrain roots (no site with high forb coverage). Make same as graminoids?
    # Leaving it alone for now.
    printnote('Not enough forb biomass at any site to estimate associated root biomass. Assume the ratio is the same as grasses?')
    printnote('No SLA measurements for forbs. Current forb value is {forbsla:1.2g}'.format(forbsla=params['slatop'].values[pft_names.index('arctic_forb')])) 
    printnote('Should forbs be evergreen or deciduous?')
    printnote('Should forbs be given quasi-woody rhizomes like grasses?')

    pft='arctic_forb' 
    change_param('leafcn',pft,obs_leafCN[:,'forb'].mean())
    change_param('frootcn',pft,obs_frootCN)
 
    change_param('fcur',pft,fcur_deciduous)
    change_param('froot_leaf',pft,1.5)
    change_param('flnr',pft,0.2)
        
    change_param('stem_leaf',pft,0.1)
    change_param('livewdcn',pft,50)   # If this is zero the model crashes on initialization
    change_param('rhizome_long',pft,rhizome_lifetime)

    # Lichen
    change_param('froot_leaf','arctic_lichen',0.2)
    change_param('leaf_long','arctic_lichen',10.0)
    change_param('roota_par','arctic_lichen',400.0)
    change_param('lflitcn','arctic_lichen',115.0)
    change_param('leafcn','arctic_lichen',84.0)
    change_param('slatop','arctic_lichen',0.036)

    # Bryophyte
    change_param('leaf_long','arctic_bryophyte',5.0)
    change_param('roota_par','arctic_bryophyte',100.0)
    change_param('lflitcn','arctic_bryophyte',66.0)
    change_param('leafcn','arctic_bryophyte',55.0)
    change_param('slatop','arctic_bryophyte',0.061)

    change_param('lflitcn','arctic_lichen',   115.0)
    change_param('lflitcn','arctic_bryophyte' , 66.0)
    change_param('lflitcn','arctic_evergreen_shrub_dwarf'    , 56.0)
    change_param('lflitcn','arctic_evergreen_shrub_tall'     , 56.0)
    change_param('lflitcn','arctic_deciduous_shrub_dwarf'    , 71.0)
    change_param('lflitcn','arctic_deciduous_shrub_low'      , 71.0)
    change_param('lflitcn','arctic_deciduous_shrub_tall'     , 71.0)
    change_param('lflitcn','arctic_deciduous_shrub_alder'    , 71.0)
    change_param('lflitcn','arctic_forb'                     , 31.0)
    change_param('lflitcn','arctic_dry_graminoid'            , 48.0)
    change_param('lflitcn','arctic_wet_graminoid'            , 65.0)

    change_param('rootb_par','arctic_lichen',   800.0)
    change_param('rootb_par','arctic_bryophyte'             ,   200.0)
    change_param('rootb_par','arctic_evergreen_shrub_dwarf' ,    13.0)
    change_param('rootb_par','arctic_evergreen_shrub_tall'  ,    13.0)
    change_param('rootb_par','arctic_deciduous_shrub_dwarf' ,    10.0)
    change_param('rootb_par','arctic_deciduous_shrub_low'   ,    10.0)
    change_param('rootb_par','arctic_deciduous_shrub_tall'  ,    10.0)
    change_param('rootb_par','arctic_deciduous_shrub_alder' ,    10.0)
    change_param('rootb_par','arctic_forb'                  ,     9.0)
    change_param('rootb_par','arctic_dry_graminoid'         ,     9.0)
    change_param('rootb_par','arctic_wet_graminoid'         ,     9.0)

    change_param('roota_par','arctic_lichen',   400.0)
    change_param('roota_par','arctic_bryophyte'             ,   200.0)
    change_param('roota_par','arctic_evergreen_shrub_dwarf' ,    30.0)
    change_param('roota_par','arctic_evergreen_shrub_tall'  ,    30.0)
    change_param('roota_par','arctic_deciduous_shrub_dwarf' ,    13.0)
    change_param('roota_par','arctic_deciduous_shrub_low'   ,    13.0)
    change_param('roota_par','arctic_deciduous_shrub_tall'  ,    13.0)
    change_param('roota_par','arctic_deciduous_shrub_alder' ,    13.0)
    change_param('roota_par','arctic_forb'                  ,     11.0)
    change_param('roota_par','arctic_dry_graminoid'         ,     11.0)
    change_param('roota_par','arctic_wet_graminoid'         ,     11.0)

    # Set up params for new dormant maintenance respiration
    # Verity suggests temperature threshold should be more like -2.5, based on Monson et al 2006, instead of +2.5
    dormant_mr_temp=273.15-1.0
    dormant_mr_factor=5e-2
    printnote('Setting dormancy temperature to {0:1.1f} C'.format(dormant_mr_temp-273.15))
    printnote('Setting dormancy maintenance resp factor to {0:1.1g}'.format(dormant_mr_factor))
    params['dormant_mr_temp']=xarray.DataArray(name='dormant_mr_temp',dims='allpfts',data=[dormant_mr_temp],attrs={'units':'degrees K','long_name':'Maximum temperature for dormant maintenance respiration'})
    params['dormant_mr_factor']=xarray.DataArray(name='dormant_mr_factor',dims='allpfts',data=[dormant_mr_factor],attrs={'units':'unitless','long_name':'Dormant maintenance respiration multiplication factor'})

    # Set up N fixation params
    params['Nfix_NPP_c1']=xarray.DataArray(name='Nfix_NPP_c1',dims='pft',data=zeros(len(pft_names))+1.8,attrs={'units':'gN/m2/yr','long_name':'Pre-exponential factor in NPP-Nfix equation'})
    params['Nfix_NPP_c2']=xarray.DataArray(name='Nfix_NPP_c2',dims='pft',data=zeros(len(pft_names))+0.003,attrs={'units':'gN/m2/yr','long_name':'Exponential factor in NPP-Nfix equation'})
    change_param('Nfix_NPP_c1','arctic_deciduous_shrub_alder',10.0)
    change_param('Nfix_NPP_c1','arctic_forb',5.0)
    change_param('Nfix_NPP_c2','arctic_deciduous_shrub_alder',0.05)

    # After checking the code, this appears to be ignored in favor of an equation based on annual average temperature
    change_universal_param('crit_onset_gdd',250)

    print('Saving params file to clm_params_updated.nc')
    params.to_netcdf(basedir+'/clm_params_updated.nc',format='NETCDF4_CLASSIC') # netcdf4 engine hangs on cades when netcdf4 version >1.3.1
    



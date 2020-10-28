from pylab import *
import xarray
import numpy
import pandas

basedir='..'

params_file=basedir+'/param_files/clm_params_newpfts_c180524_orig.nc'
params_file_default=basedir+'/param_files/clm_params_defaultpfts_c180524_orig.nc'

params_fengming=xarray.open_dataset(params_file)
params_default=xarray.open_dataset(params_file_default)
params_new=xarray.open_dataset(basedir+'/clm_params_updated.nc')

pft_names=[name.strip() for name in params_fengming['pftname'].values.astype(str)]
pft_names_default=[name.strip() for name in params_default['pftname'].values.astype(str)]

def prettify_pft_name(name):
    if name.startswith('arctic'):
        pretty_name=name[len('arctic_'):]
    else:
        pretty_name=name
    pretty_name = ' '.join(pretty_name.split('_')).title()
    if pretty_name == "Dry Graminoid":
        pretty_name = "Graminoid"
    if len(pretty_name.split())>=2 and pretty_name.split()[-2]=='Shrub':
        l=pretty_name.split()
        end=l.pop(-1)
        l.insert(-1,end)
        pretty_name=' '.join(l)
    if pretty_name == "Deciduous Tall Shrub":
        pretty_name="Deciduous Low to Tall Shrub"
    return pretty_name

def get_time(data):
    if isinstance(columndata['time'].data[0],numpy.datetime64):
        t=pandas.DatetimeIndex(columndata['time'].data)
        return t.year+(t.dayofyear-1)/365
    else:
        return array([tt.year + (tt.dayofyr-1)/365 for tt in columndata['time'].data])

# domaindata=xarray.open_dataset('/lustre/or-hydra/cades-ccsi/proj-shared/project_acme/cesminput_ngee/share/domains/domain.clm/domain.lnd.51x63pt_kougarok-NGEE_TransA_navy.nc')
# surfdata=xarray.open_dataset(basedir+'/param_files/surfdata_51x63pt_kougarok-NGEE_TransA_simyr1850_c181115-sub12.nc')
# surfdata=xarray.open_dataset(basedir+'/param_files/surfdata_51x63pt_kougarok-NGEE_TransA_simyr1850_c181115-sub12_updated_2019-02-15.nc')
# surfdata=xarray.open_dataset(basedir+'/param_files/surfdata_51x63pt_kougarok-NGEE_TransA_simyr1850_c190604-sub12_updated_2019-06-17.nc')
surfdata=xarray.open_dataset(basedir+'/param_files/surfdata_51x63pt_kougarok-NGEE_TransA_simyr1850_c190604-sub12_updated_2020-10-08.nc')
# surfdata_default=xarray.open_dataset(basedir+'/param_files/surfdata_51x63pt_kougarok-NGEE_TransA_simyr1850_c181115default.nc')
surfdata_default=xarray.open_dataset(basedir+'/param_files/surfdata_Kougarok_defaultPFTs_all-shrubs-decid-boreal.nc')

# landscape_ecotypes=['NAMC','DSLT','AS','WBT','TTWBT','TT']
# ecotype_names={'NAMC':'Non-acidic mountain complex',
#                'DSLT':'Dwarf shrub lichen tundra',
#                'AS'  :'Alder shrubland',
#                'WBT' :'Willow birch tundra',
#                'TTWBT':'Tussock tundra/willow birch tundra',
#                'TT'  :'Tussock tundra'}
# landscape_ecotypes=['DL','DSLT','AS','WBT','ASV','TT']
# ecotype_names={'DL':'Dryas-lichen dwarf shrub tundra',
#                'DSLT':'Dwarf shrub lichen tundra',
#                'AS'  :'Alder shrubland',
#                'WBT' :'Willow birch tundra',
#                'ASV':'Alder savanna',
#                'TT'  :'Tussock tundra'}
landscape_ecotypes=['DLST','BEL','AS','WBT','ASV','TT']
ecotype_names={'DLST':'Dryas-lichen dwarf shrub tundra',
               'BEL':'Birch-Ericaceous-lichen shrub tundra',
               'AS'  :'Alder shrubland',
               'WBT' :'Willow birch tundra',
               'ASV':'Alder savanna',
               'TT'  :'Tussock tundra',
               'Gridcell':'Downscaled grid cell'}


obsdata_PFT_mappings={'dwarf shrub deciduous':'arctic_deciduous_shrub_dwarf',
                        'deciduous dwarf shrub':'arctic_deciduous_shrub_dwarf',
                      'dwarf shrub evergreen':'arctic_evergreen_shrub_dwarf',
                      'evergreen shrub':'arctic_evergreen_shrub_dwarf',
                      'forb':'arctic_forb',
                      'graminoid':'arctic_dry_graminoid', # **** model has wet and dry graminoids
                      'lichen':'arctic_lichen',
                      'low shrub deciduous':'arctic_deciduous_shrub_low',
                      'deciduous low shrub':'arctic_deciduous_shrub_low',
                      'mixed':'not_vegetated',            # *****
                      'moss':'arctic_bryophyte',
                      'bryophyte':'arctic_bryophyte',
                        'other':'not_vegetated',          # ***** what to do with this?
                        'tall shrub deciduous alder':'arctic_deciduous_shrub_alder',
                        'tall shrub deciduous birch':'arctic_deciduous_shrub_tall', # **** Model is not separating birch and willow
                        'tall shrub deciduous willow':'arctic_deciduous_shrub_tall', # **** Model is not separating birch and willow
                        'potential tall shrub deciduous alder':'arctic_deciduous_shrub_alder',
                        'potential tall shrub deciduous birch':'arctic_deciduous_shrub_tall', # **** Model is not separating birch and willow
                        'potential tall shrub deciduous willow':'arctic_deciduous_shrub_tall', # **** Model is not separating birch and willow
                        'potential tall shrub deciduous non-alder':'arctic_deciduous_shrub_tall',
                        'deciduous low to tall willow and birch shrub':'arctic_deciduous_shrub_tall',
                        'deciduous low to tall alder shrub':'arctic_deciduous_shrub_alder',
                        'non-vegetated':'not_vegetated'
                        }
                        
obsdata_defaultPFT_mappings={'dwarf shrub deciduous':'broadleaf_deciduous_boreal_shrub',
                      'dwarf shrub evergreen':'broadleaf_evergreen_shrub',
                      'forb':'c3_arctic_grass',
                      'graminoid':'c3_arctic_grass', # **** model has wet and dry graminoids
                      'lichen':'nonvascular',
                      'low shrub deciduous':'broadleaf_deciduous_boreal_shrub',
                      'mixed':'not_vegetated',            # *****
                      'moss':'nonvascular',
                      'bryophyte':'nonvascular',
                        'other':'not_vegetated',          # ***** what to do with this?
                        'tall shrub deciduous alder':'broadleaf_deciduous_boreal_shrub',
                        'tall shrub deciduous birch':'broadleaf_deciduous_boreal_shrub', # **** Model is not separating birch and willow
                        'tall shrub deciduous willow':'broadleaf_deciduous_boreal_shrub', # **** Model is not separating birch and willow
                        'potential tall shrub deciduous alder':'broadleaf_deciduous_boreal_shrub',
                        'potential tall shrub deciduous birch':'broadleaf_deciduous_boreal_shrub', # **** Model is not separating birch and willow
                        'potential tall shrub deciduous willow':'broadleaf_deciduous_boreal_shrub', # **** Model is not separating birch and willow
                        'potential tall shrub deciduous non-alder':'broadleaf_deciduous_boreal_shrub'
                        }

obsdata_E3SMPFT_mappings={'dwarf shrub deciduous':'broadleaf_deciduous_boreal_shrub',
                      'dwarf shrub evergreen':'broadleaf_deciduous_boreal_shrub',
                      'forb':'c3_arctic_grass',
                      'graminoid':'c3_arctic_grass', # **** model has wet and dry graminoids
                      'lichen':'nonvascular',
                      'low shrub deciduous':'broadleaf_deciduous_boreal_shrub',
                      'mixed':'not_vegetated',            # *****
                      'moss':'nonvascular',
                      'bryophyte':'nonvascular',
                        'other':'not_vegetated',          # ***** what to do with this?
                        'tall shrub deciduous alder':'broadleaf_deciduous_boreal_shrub',
                        'tall shrub deciduous birch':'broadleaf_deciduous_boreal_shrub', # **** Model is not separating birch and willow
                        'tall shrub deciduous willow':'broadleaf_deciduous_boreal_shrub', # **** Model is not separating birch and willow
                        'potential tall shrub deciduous alder':'broadleaf_deciduous_boreal_shrub',
                        'potential tall shrub deciduous birch':'broadleaf_deciduous_boreal_shrub', # **** Model is not separating birch and willow
                        'potential tall shrub deciduous willow':'broadleaf_deciduous_boreal_shrub', # **** Model is not separating birch and willow
                        'potential tall shrub deciduous non-alder':'broadleaf_deciduous_boreal_shrub'
                        }

default_new_mappings={'broadleaf_evergreen_shrub':'arctic_evergreen_shrub_dwarf',
                      'broadleaf_deciduous_boreal_shrub':'arctic_deciduous_shrub_tall',
                      'c3_arctic_grass':'arctic_dry_graminoid'}



pft_colors_dict={
 'not_vegetated':'0.9',
 'arctic_lichen':'#FF8236',
 'arctic_bryophyte':'#BD2D2D',
 'arctic_evergreen_shrub_dwarf':'#FFE802',
 'arctic_evergreen_shrub_tall':'w',
 'arctic_deciduous_shrub_dwarf':'#FF4082',
 'arctic_deciduous_shrub_low':'#46C7C7',
 'arctic_deciduous_shrub_tall':'#4D9E6E',
 'arctic_deciduous_shrub_alder':'#354A24',
 'arctic_forb':'#BB80FF',
 'arctic_dry_graminoid':'#44E300',
 'arctic_wet_graminoid':'#44E300',
 
}
pft_colors=[pft_colors_dict[name] for name in pft_names]
pft_colors_default=['0.9' for n in range(len(pft_names_default))]
pft_colors_default[pft_names_default.index('broadleaf_evergreen_shrub')]='C3'
pft_colors_default[pft_names_default.index('broadleaf_deciduous_boreal_shrub')]='0.3'
pft_colors_default[pft_names_default.index('c3_arctic_grass')]=pft_colors_dict['arctic_dry_graminoid']
pft_colors_default[pft_names_default.index('needleleaf_evergreen_boreal_tree')]='darkgreen'

PFT_percents=pandas.DataFrame(data=surfdata.PCT_NAT_PFT.values.squeeze(),index=pft_names,columns=landscape_ecotypes)
PFT_percents_default=pandas.DataFrame(data=surfdata_default.PCT_NAT_PFT.values.squeeze(),index=pft_names_default[:17],columns=landscape_ecotypes)

# Read directly from Amy's spreadsheet to make sure model relative areas are the same as her latest measurements
PFT_percents_latest=pandas.read_excel('../obs_data/ngee_arctic_pft_tallies_kougarok_20191010.xlsx',sheet_name='KG_bar_graphs',header=3,nrows=10,usecols=range(7),index_col=0).rename(
    columns={
    'Dryas-lichen dwarf shrub tundra':'DLST',
    'Birch-Ericaceous lichen shrub tundra':'BEL',
    'Alder shrubland':'AS',
    'Willow-birch tundra':'WBT',
    'Alder savanna':'ASV',
    'Tussock-lichen tundra':'TT'
},index=obsdata_PFT_mappings)
PFT_percents_latest.loc['arctic_evergreen_shrub_tall']=0 
PFT_percents_latest.loc['arctic_wet_graminoid']=0 
PFT_percents_latest['ASV']['arctic_lichen']=PFT_percents_latest['ASV']['arctic_lichen']+0.01
PFT_percents_latest=PFT_percents_latest.loc[pft_names]

surfdata_E3SM=xarray.open_dataset('../param_files/surfdata_Kougarok_downscaled_PFTs_soildepths.nc')
PFT_percents_E3SM=pandas.DataFrame(data=surfdata_E3SM.PCT_NAT_PFT.values.squeeze(),index=pft_names_default[:17],columns=landscape_ecotypes)

def read_pftfile(filename,maxyear=None,decode_times=True):
    with xarray.open_dataset(filename,decode_times=decode_times) as output_PFTs:
        if maxyear is not None:
            output_PFTs=output_PFTs.sel(time=array([xx.year for xx in output_PFTs.time.values])<=maxyear)

        pft_mask=output_PFTs.pfts1d_itype_lunit == 1
        weights=output_PFTs.pfts1d_wtgcell[pft_mask]
        if len(weights) == len(pft_names_default[:17])*len(output_PFTs.lndgrid):
            print('Reading from sim with default PFTs')
            vegdata_PFTs=xarray.Dataset(coords={'time':output_PFTs.time,'PFT':arange(len(pft_names_default[:17])),'ecotype':arange(len(output_PFTs.lndgrid))})
            newshape=(len(output_PFTs.lndgrid),len(pft_names_default[:17]))
            vegdata_PFTs.coords['PFTnames']=xarray.DataArray(pft_names_default[:17],dims=('PFT',)) 
            vegdata_PFTs.coords['PFTcolors']=xarray.DataArray(pft_colors_default[:17],dims=('PFT',))
        else:
            print('Reading from sim with new PFTs')
            vegdata_PFTs=xarray.Dataset(coords={'time':output_PFTs.time,'PFT':arange(len(pft_names)),'ecotype':arange(len(output_PFTs.lndgrid))})
            newshape=(len(output_PFTs.lndgrid),len(pft_names))
            vegdata_PFTs.coords['PFTnames']=xarray.DataArray(pft_names,dims=('PFT',))
            vegdata_PFTs.coords['PFTcolors']=xarray.DataArray(pft_colors,dims=('PFT',))

        vegdata_PFTs['weights']=(('ecotype','PFT'),weights.values.reshape(newshape))
        if len(output_PFTs.lndgrid)>1:
            vegdata_PFTs.coords['community_names']=xarray.DataArray(landscape_ecotypes,dims=('ecotype',))
        else:
            vegdata_PFTs.coords['community_names']=xarray.DataArray(['Gridcell'],dims=('ecotype',))


        for varname in output_PFTs.variables:
            var=output_PFTs[varname]
            if var.dims == ('time','pft'):
                print(varname)
                vardata=var.values[:,pft_mask].reshape((len(var.time),newshape[0],newshape[1]))
                vegdata_PFTs[var.name+'_unweighted']=(('time','ecotype','PFT'),vardata)
                vegdata_PFTs[var.name+'_unweighted'].attrs['long_name']=var.long_name
                vegdata_PFTs[var.name+'_unweighted'].attrs['units']=var.units
            elif var.dims == ('time','levdcmp','pft'):
                print(varname+' (depth-resolved)')
                vardata=var.values[:,:,pft_mask].reshape((len(var.time),len(var.levdcmp),newshape[0],newshape[1]))
                vegdata_PFTs[var.name+'_unweighted']=(('time','levdcmp','ecotype','PFT'),vardata)
                vegdata_PFTs[var.name+'_unweighted'].attrs['long_name']=var.long_name
                vegdata_PFTs[var.name+'_unweighted'].attrs['units']=var.units
            else:
                print('NonPFT variable %s'%varname)
                vegdata_PFTs[var.name]=var
        return vegdata_PFTs


ecotype_names_list=[ecotype_names[landscape_ecotypes[n]] for n in surfdata.lsmlon.values]


def plot_var(varname,obsdata=None,minyear=2008,maxyear=2010,longname=None,units=None):
    if isinstance(varname,str):
        dat_default=data_default[varname]
        dat=data[varname]
        if longname is None:
            longname=data[varname].long_name
        if units is None:
            units=data[varname].units
    else:
        dat_default=data_default[varname[0]]
        dat=data[varname[0]]
        if len(varname)>1:
            if longname is None:
                raise ValueError('Must supply longname if passing multiple data fields')
            for name in varname[1:]:
                dat_default=dat_default+data_default[name]
                dat=dat+data[name]
        if longname is None:
            longname=data[varname[0]].long_name
        if units is None:
            units=data[varname[0]].units


    import datetime
    mindate=datetime.date(minyear,1,1)
    maxdate=datetime.date(maxyear,1,1)
    defaultparams_ax=subplot(211)
    title('Default PFTs')
    plot(t_default,dat_default)
    xlabel('Time')
    ylabel('%s (%s)'%(longname,units))
    defaultparams_ax.set_xlim(left=mindate,right=maxdate)


    newparams_ax=subplot(212)
    title('Kougarok PFTs')
    plot(t,dat)
    legend(labels=ecotype_names_list,ncol=2,fontsize='small')
    ylabel('%s (%s)'%(longname,units))
    xlabel('Time')
    newparams_ax.set_xlim(left=datetime.date(minyear,1,1),right=datetime.date(maxyear,1,1))

    if obsdata is not None:
        # Collection was in 2016-2017, but model only goes through 2010. No big deal I guess.
        for ecotype in landscape_ecotypes:
            defaultparams_ax.plot([mindate,maxdate],[obsdata[ecotype],obsdata[ecotype]],c='C%d'%(landscape_ecotypes.index(ecotype)),ls='--')
            newparams_ax.plot([mindate,maxdate],[obsdata[ecotype],obsdata[ecotype]],c='C%d'%(landscape_ecotypes.index(ecotype)),ls='--')

    tight_layout()
    return defaultparams_ax,newparams_ax

# Older dataset
Koug_meas_biomass_old=pandas.read_excel(basedir+'/obs_data/NGEEArctic_Q3ELM_KougarokBiomass&NPP_20181112.xlsx',sheet_name='data')\
    .set_index(['Ecotype','ELMgroup'])
Koug_meas_chem_old=pandas.read_excel(basedir+'/obs_data/NGEEArctic_Q3ELM_KougarokSLA&Chemistry_20181112.xlsx',sheet_name='data')\
    .rename(columns={'ELM_PFT':'ELMgroup'}).set_index(['Ecotype','ELMgroup'])

# New version
# Koug_meas_biomass=pandas.read_excel(basedir+'/obs_data/NGEE Arctic Veg data compiled 12June2019/NGEEArctic_Q3ELM_KougarokBiomass&NPP_20190612.xlsx',sheet_name='data')\
#     .set_index(['Ecotype','PlotID','ELM_PFT'])
# Koug_meas_chem=pandas.read_excel(basedir+'/obs_data/NGEE Arctic Veg data compiled 12June2019/Kougarok_Q3ELM_KougarokSLA&Chemistry_20190612.xlsx',sheet_name='data')\
#     .set_index(['Ecotype','PlotID','ELM_PFT'])

Koug_meas_biomass=pandas.read_excel(basedir+'/obs_data/NGEE Arctic Veg data compiled 4May2020/NGEEArctic_Q3ELM_Biomass&NPP_20200504.xlsx',sheet_name='final data')\
    .set_index(['Ecotype','PlotID','ELM_PFT'])
Koug_meas_chem=pandas.read_excel(basedir+'/obs_data/NGEE Arctic Veg data compiled 4May2020/Kougarok_Q3ELM_KougarokSLA&Chemistry_20200504.xlsx',sheet_name='data')\
    .set_index(['Ecotype','PlotID','ELM_PFT'])

# Data from Amy measurements: doi:10.5440/1465967
Breen_data=pandas.read_csv('../obs_data/ngee_arctic_kougarok_2016_veg_comp_env_table_v1_20180828.csv',header=0,skiprows=[1,2],skipinitialspace=True).dropna(axis='index',how='all')
obs_heights=pandas.DataFrame(index=landscape_ecotypes)
ecotype_height_mapping={
    'Moist to dry alder (Alnus viridis) communities and alder savannas':'AS',
    'Zonal habitats with erect-dwarf-shrub tundra acidic soils, subzones D and E':'BEL',
    'Dry azonal habitats, base-rich soils, subzones D and E':'DLST',
    'Tussock tundra (Eriophorum vaginatum)':'TT',
    'Moist to dry alder (Alnus viridis) communities and alder savannas':'AS',
    'Low-shrub tundra, acidic soils, warmest parts of subzone E':'WBT'
}
ecotype_height_mapping={
    'Alder shrubland':'AS',
    'Dwarf-shrub lichen tundra':'BEL',
    'Non-acidic mountain complex':'DLST',
    'Tussock tundra':'TT',
    'Tussock tundra-Willow-birch tundra complex OR Alder savannah in tussock tundra':'ASV',
    'Willow-birch tundra':'WBT'
}
Breen_means=Breen_data.groupby('preliminary_plant_community_name').mean().rename(index=ecotype_height_mapping)  
obs_heights['tree_height_mean']=Breen_means['mean_tree_layer_height']
obs_heights['shrub_height_mean']=Breen_means['mean_shrub_layer_height']
obs_heights['tall_shrub_height_mean']=Breen_means['mean_tall_shrub_height']
obs_heights['low_shrub_height_mean']=Breen_means['mean_low_shrub_height']
obs_heights['dwarf_shrub_height_mean']=Breen_means['mean_dwarf_shrub_height']
obs_heights['forb_height_mean']=Breen_means['mean_forb_height']

obs_heights['canopy_height_max']=Breen_data.groupby('preliminary_plant_community_name').max().rename(index=ecotype_height_mapping)['maximum_ canopy_height']

soildepth=pandas.read_excel('../obs_data/NGEEArctic_Kougarok_SoilDepth_starting2016_v1.xlsx',header=6)
soilcores=pandas.read_excel('../obs_data/NGEE Arctic_Kougarok_2016_Soil cores_20170906.xlsx')

# convert from cm to m
obs_heights=obs_heights/100

def plot_var_PFTs(varname,moddata,ecotype_num,ax=None,obsdata=None,minyear=0,maxyear=150,longname=None,units=None,modfactor=1.0,cumulative=False,weight_area=True,plotsum=False,**kwargs):
    if isinstance(varname,str):
        dat=moddata[varname+'_unweighted'].copy()
        if longname is None:
            longname=dat.long_name
        if units is None:
            units=dat.units
    else:
        # Assuming we are adding multiple data fields together
        dat=moddata[varname[0]+'_unweighted'].copy()
        if len(varname)>1:
            if longname is None:
                raise ValueError('Must supply longname if passing multiple data fields')
            units=dat.units
            for name in varname[1:]:
                dat=dat+moddata[name+'_unweighted']
        else:
            if longname is None:
                longname=dat.long_name
            if units is None:
                units=dat.units

    if weight_area:
        dat=dat*moddata.weights

    import datetime

    if isinstance(moddata['time'].data[0],numpy.number):
        t=moddata['time']/365
        dt=diff(t).mean()
        mindate=minyear
        maxdate=maxyear
        if cumulative:
            dat=dat.cumsum(dim='time')*dt
            dat=dat-dat.isel(time=nonzero(t>minyear)[0][0])
            
    elif not isinstance(moddata['time'].data[0],numpy.datetime64):
        t=array([tt.year + (tt.dayofyr-1)/365 for tt in moddata['time'].data])
        dt=diff(t).mean()

        mindate=minyear
        maxdate=maxyear
        if cumulative:
            dat=dat.cumsum(dim='time')*dt
            dat=dat-dat.isel(time=nonzero(t>minyear)[0][0])

    else:
        t=moddata['time'].data
        dt=diff(t).mean()
        mindate=datetime.date(minyear,1,1)
        maxdate=datetime.date(maxyear,1,1)

        if cumulative:
            dat=dat.cumsum(dim='time')*dt
            dat=dat-dat.sel(time=mindate)


    if ax is None:
        ax=gca()

    dat=dat*modfactor
    ax.set_title(longname)
    for pft in moddata.PFT.values:
        pftlabel=str(moddata.PFTnames.sel(PFT=pft).values)
        c=str(moddata.PFTcolors.sel(PFT=pft).values)
        if pftlabel.startswith('arctic_'):
            pftlabel=pftlabel[len('arctic_'):]
        ax.plot(t,dat.sel(ecotype=ecotype_num,PFT=pft),c=c,label=pftlabel,**kwargs)
    if plotsum:
        ax.plot(t,dat.sel(ecotype=ecotype_num).sum(dim='PFT'),c='C0',label='Sum of PFTs',**kwargs)
    ax.set_xlabel('Time')
    ax.set_ylabel('%s (%s)'%(longname,units))
    ax.set_xlim(left=mindate,right=maxdate)
    ax.xaxis.set_tick_params(rotation=0)

    if obsdata is not None:
        # Collection was in 2016-2017, but model only goes through 2010. No big deal I guess.
        if 'ELM_PFT' in obsdata.index.names:
            if 'PlotID' in obsdata.index.names:
                obsdata_mean=obsdata[landscape_ecotypes[ecotype_num]].mean(level='ELM_PFT')
                obsdata_std=obsdata[landscape_ecotypes[ecotype_num]].std(level='ELM_PFT')
            else:
                obsdata_mean=obsdata[landscape_ecotypes[ecotype_num]]
                obsdata_std=nan
            for pft in obsdata_mean.index:
                ax.fill_between([mindate,maxdate],[(obsdata_mean-obsdata_std)[pft],(obsdata_mean-obsdata_std)[pft]],
                            [(obsdata_mean+obsdata_std)[pft],(obsdata_mean+obsdata_std)[pft]],
                            color=pft_colors[(pft_names.index(obsdata_PFT_mappings[pft]))],alpha=0.3)
                # ax.plot([mindate,maxdate],[obsdata_mean,obsdata_mean],c=pft_colors[(pft_names.index(obsdata_PFT_mappings[pft]))],ls='--')
        else:
            if 'PlotID' in obsdata.index.names:
                obsdata_mean=obsdata[landscape_ecotypes[ecotype_num]].mean()
                obsdata_std=obsdata[landscape_ecotypes[ecotype_num]].std()
            else:
                obsdata_mean=obsdata[landscape_ecotypes[ecotype_num]]
                obsdata_std=nan

            ax.fill_between([mindate,maxdate],[(obsdata_mean-obsdata_std),(obsdata_mean-obsdata_std)],
                        [(obsdata_mean+obsdata_std),(obsdata_mean+obsdata_std)],
                        color=pft_colors[pft_names.index('not_vegetated')],alpha=0.3)
            # ax.plot([mindate,maxdate],[obsdata_mean,obsdata_mean],c=pft_colors[pft_names.index('not_vegetated')],ls='--')

    return dat


def plot_pair(var,vegdata_oldparams,vegdata_newparams,ecotype_num=1,**kwargs):
    fig=figure(var,figsize=(6.4,6.4))
    fig.clf()
    nplots=2
    gs=fig.add_gridspec(ncols=1,nrows=2)
    subplot_handles={}
    subplot_handles[var+'_old']=fig.add_subplot(gs[0,0])
    subplot_handles[var+'_new']=fig.add_subplot(gs[1,0])

    plot_var_PFTs(var,vegdata_oldparams,ecotype_num,ax=subplot_handles[var+'_old'],**kwargs)
    plot_var_PFTs(var,vegdata_newparams,ecotype_num,ax=subplot_handles[var+'_new'],**kwargs)
    leg=subplot_handles[var+'_old'].legend(fontsize='small')
    leg.set_draggable(True)

    tight_layout()
    return fig
    
def get_var_PFTs(varname,moddata,weight_area=True):
    if isinstance(varname,str):
        dat=moddata[varname+'_unweighted'].copy()
    else:
        # Assuming we are adding multiple data fields together
        dat=moddata[varname[0]+'_unweighted'].copy()
        if len(varname)>1:
            for name in varname[1:]:
                dat=dat+moddata[name+'_unweighted']


    if weight_area:
        dat=dat*moddata.weights
        
    return dat


def save_all_figs(dirname=basedir+'/Figures',format='png',**kwargs):
    for fname in get_figlabels():
        fname_fixed=fname.replace('/','-')
        print(fname_fixed)
        figure(fname_fixed).savefig('{dirname:s}/{fname:s}.{format}'.format(dirname=dirname,format=format,fname=fname_fixed),**kwargs)


def pft_params(paramdata,paramnames):
    pftnames=[name.strip() for name in paramdata['pftname'].values.astype(str)]
    if isinstance(paramnames,str):
        paramnames=[paramnames]
    pdict={'pftname':pftnames}
    for name in paramnames:
        pdict[name]=paramdata[name]
    return pandas.DataFrame(pdict).set_index('pftname')
    
def plot_PFT_distributions(axs=None):
    if axs is None:
        axs=gcf().subplots(ncols=3)
    sca(axs[0])
    names=[]
    bottom=zeros(len(landscape_ecotypes))
    for pftnum in range(len(pft_names_default[:17])):
        pft_pcts=PFT_percents_E3SM.loc[pft_names_default[pftnum]]
        if (pft_pcts==0).all():
            continue
        bar(arange(len(landscape_ecotypes)),pft_pcts,bottom=bottom,facecolor=pft_colors_default[pftnum])
        bottom=bottom+pft_pcts
        names.append(prettify_pft_name(pft_names_default[pftnum] ))

    xticks(arange(len(landscape_ecotypes)),landscape_ecotypes,rotation=0)
    title('E3SM grid PFTs')
    l=legend(labels=names,loc=(0.0,1.1),fontsize='small')
    l.set_draggable(True)
    
    sca(axs[1])
    names=[]
    bottom=zeros(len(landscape_ecotypes))
    for pftnum in range(len(pft_names_default[:17])):
        pft_pcts=PFT_percents_default.loc[pft_names_default[pftnum]]
        if (pft_pcts==0).all():
            continue
        bar(arange(len(landscape_ecotypes)),pft_pcts,bottom=bottom,facecolor=pft_colors_default[pftnum])
        bottom=bottom+pft_pcts
        names.append(prettify_pft_name(pft_names_default[pftnum] ))

    xticks(arange(len(landscape_ecotypes)),landscape_ecotypes,rotation=0)
    title('Default ELM PFTs, site areas')
    l=legend(labels=names,loc=(0.0,1.1),fontsize='small')
    l.set_draggable(True)

    sca(axs[2])
    names=[]
    bottom=zeros(len(landscape_ecotypes))
    order=['not_vegetated','arctic_deciduous_shrub_dwarf','arctic_evergreen_shrub_dwarf','arctic_deciduous_shrub_low',
            'arctic_deciduous_shrub_tall','arctic_deciduous_shrub_alder',
            'arctic_forb','arctic_dry_graminoid','arctic_bryophyte','arctic_lichen']
    # for pftnum in range(len(pft_names)):
    for pftname in order:
        pftnum=pft_names.index(pftname)
        pft_pcts=PFT_percents.loc[pft_names[pftnum]]
        if (pft_pcts==0).all():
            continue
        bar(arange(len(landscape_ecotypes)),pft_pcts,bottom=bottom,facecolor=pft_colors[pftnum])
        bottom=bottom+pft_pcts
        names.append(prettify_pft_name(pft_names[pftnum] ))
        
    # Add indicator of shrubs and non-shrubs
    notveg=PFT_percents.loc['not_vegetated']
    allshrubs=PFT_percents[PFT_percents.index.str.contains('shrub')].sum()
    nonshrubs=100-notveg-allshrubs
    bar(arange(len(landscape_ecotypes)),notveg,facecolor='None',edgecolor=pft_colors_default[pft_names_default.index('not_vegetated')],linewidth=2,linestyle='--',alpha=0.8)
    bar(arange(len(landscape_ecotypes)),allshrubs,bottom=notveg,facecolor='None',edgecolor=pft_colors_default[pft_names_default.index('broadleaf_deciduous_boreal_shrub')],linewidth=1.5,linestyle='--',alpha=0.8)
    bar(arange(len(landscape_ecotypes)),nonshrubs,bottom=notveg+allshrubs,facecolor='None',edgecolor=pft_colors_default[pft_names_default.index('c3_arctic_grass')],linewidth=1.5,linestyle='--',alpha=0.8)

    xticks(arange(len(landscape_ecotypes)),landscape_ecotypes,rotation=0)
    title('Arctic PFTs')
    l=legend(labels=names,loc=(0.0,1.1),fontsize='small',ncol=2)
    l.set_draggable(True)

    # tight_layout()

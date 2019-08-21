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



# domaindata=xarray.open_dataset('/lustre/or-hydra/cades-ccsi/proj-shared/project_acme/cesminput_ngee/share/domains/domain.clm/domain.lnd.51x63pt_kougarok-NGEE_TransA_navy.nc')
# surfdata=xarray.open_dataset(basedir+'/param_files/surfdata_51x63pt_kougarok-NGEE_TransA_simyr1850_c181115-sub12.nc')
surfdata=xarray.open_dataset(basedir+'/param_files/surfdata_51x63pt_kougarok-NGEE_TransA_simyr1850_c181115-sub12_updated_2019-02-15.nc')
surfdata_default=xarray.open_dataset(basedir+'/param_files/surfdata_51x63pt_kougarok-NGEE_TransA_simyr1850_c181115default.nc')

# landscape_ecotypes=['NAMC','DSLT','AS','WBT','TTWBT','TT']
# ecotype_names={'NAMC':'Non-acidic mountain complex',
#                'DSLT':'Dwarf shrub lichen tundra',
#                'AS'  :'Alder shrubland',
#                'WBT' :'Willow birch tundra',
#                'TTWBT':'Tussock tundra/willow birch tundra',
#                'TT'  :'Tussock tundra'}
landscape_ecotypes=['DL','DSLT','AS','WBT','ASV','TT']
ecotype_names={'DL':'Dryas-lichen dwarf shrub tundra',
               'DSLT':'Dwarf shrub lichen tundra',
               'AS'  :'Alder shrubland',
               'WBT' :'Willow birch tundra',
               'ASV':'Alder savanna',
               'TT'  :'Tussock tundra'}

pft_colors=['C%d'%n for n in range(10)] + ['k','purple']
pft_colors_default=['C0' for n in range(len(pft_names_default))]
pft_colors_default[pft_names_default.index('broadleaf_evergreen_shrub')]='C3'
pft_colors_default[pft_names_default.index('broadleaf_deciduous_boreal_shrub')]='C7'
pft_colors_default[pft_names_default.index('c3_arctic_grass')]='k'

PFT_percents=pandas.DataFrame(data=surfdata.PCT_NAT_PFT.values.squeeze(),index=pft_names,columns=landscape_ecotypes)
PFT_percents_default=pandas.DataFrame(data=surfdata_default.PCT_NAT_PFT.values.squeeze(),index=pft_names_default[:17],columns=landscape_ecotypes)

def read_pftfile(filename,maxyear=None):
    output_PFTs=xarray.open_dataset(filename)
    if maxyear is not None:
        output_PFTs=output_PFTs.sel(time=array([xx.year for xx in output_PFTs.time.values])<=maxyear)

    pft_mask=output_PFTs.pfts1d_itype_lunit == 1
    weights=output_PFTs.pfts1d_wtgcell[pft_mask]
    if len(weights) == len(pft_names_default[:17])*len(landscape_ecotypes):
        print('Reading from sim with default PFTs')
        vegdata_PFTs=xarray.Dataset(coords={'time':output_PFTs.time,'PFT':arange(len(pft_names_default[:17])),'ecotype':arange(len(landscape_ecotypes))})
        newshape=(len(landscape_ecotypes),len(pft_names_default[:17]))
        vegdata_PFTs.coords['PFTnames']=xarray.DataArray(pft_names_default[:17],dims=('PFT',)) 
        vegdata_PFTs.coords['PFTcolors']=xarray.DataArray(pft_colors_default[:17],dims=('PFT',))
    else:
        print('Reading from sim with new PFTs')
        vegdata_PFTs=xarray.Dataset(coords={'time':output_PFTs.time,'PFT':arange(len(pft_names)),'ecotype':arange(len(landscape_ecotypes))})
        newshape=(len(landscape_ecotypes),len(pft_names))
        vegdata_PFTs.coords['PFTnames']=xarray.DataArray(pft_names,dims=('PFT',))
        vegdata_PFTs.coords['PFTcolors']=xarray.DataArray(pft_colors,dims=('PFT',))

    vegdata_PFTs['weights']=(('ecotype','PFT'),weights.values.reshape(newshape))
    vegdata_PFTs.coords['community_names']=xarray.DataArray(landscape_ecotypes,dims=('ecotype',))
    for varname in output_PFTs.variables:
        var=output_PFTs[varname]
        if var.dims != ('time','pft'):
            continue
        else:
            print(varname)
            vardata=var.values[:,pft_mask].reshape((len(var.time),newshape[0],newshape[1]))
            vegdata_PFTs[var.name+'_unweighted']=(('time','ecotype','PFT'),vardata)
            vegdata_PFTs[var.name+'_unweighted'].attrs['long_name']=var.long_name
            vegdata_PFTs[var.name+'_unweighted'].attrs['units']=var.units
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

Koug_meas_biomass=pandas.read_excel(basedir+'/obs_data/NGEEArctic_Q3ELM_KougarokBiomass&NPP_20181112.xlsx',sheet_name='data')\
    .set_index(['Ecotype','ELMgroup'])
Koug_meas_chem=pandas.read_excel(basedir+'/obs_data/NGEEArctic_Q3ELM_KougarokSLA&Chemistry_20181112.xlsx',sheet_name='data')\
    .rename(columns={'ELM_PFT':'ELMgroup'}).set_index(['Ecotype','ELMgroup'])

obsdata_PFT_mappings={'dwarf shrub deciduous':'arctic_deciduous_shrub_dwarf',
                      'dwarf shrub evergreen':'arctic_evergreen_shrub_dwarf',
                      'forb':'arctic_forb',
                      'graminoid':'arctic_dry_graminoid', # **** model has wet and dry graminoids
                      'lichen':'arctic_lichen',
                      'low shrub deciduous':'arctic_deciduous_shrub_low',
                      'mixed':'not_vegetated',            # *****
                      'moss':'arctic_bryophyte',
                        'other':'not_vegetated',          # ***** what to do with this?
                        'tall shrub deciduous alder':'arctic_deciduous_shrub_alder',
                        'tall shrub deciduous birch':'arctic_deciduous_shrub_tall', # **** Model is not separating birch and willow
                        'tall shrub deciduous willow':'arctic_deciduous_shrub_tall' # **** Model is not separating birch and willow
                        }

default_new_mappings={'broadleaf_evergreen_shrub':'arctic_evergreen_shrub_dwarf',
                      'broadleaf_deciduous_boreal_shrub':'arctic_deciduous_shrub_tall',
                      'c3_arctic_grass':'arctic_dry_graminoid'}

# Data from Amy measurements: doi:10.5440/1465967
Breen_data=pandas.read_csv('../obs_data/ngee_arctic_kougarok_2016_veg_comp_env_table_v1_20180828.csv',header=0,skiprows=[1,2],skipinitialspace=True).dropna(axis='index',how='all')
obs_heights=pandas.DataFrame(index=landscape_ecotypes)
ecotype_height_mapping={
    'Moist to dry alder (Alnus viridis) communities and alder savannas':'AS',
    'Zonal habitats with erect-dwarf-shrub tundra acidic soils, subzones D and E':'DSLT',
    'Dry azonal habitats, base-rich soils, subzones D and E':'DL',
    'Tussock tundra (Eriophorum vaginatum)':'TT',
    'Moist to dry alder (Alnus viridis) communities and alder savannas':'AS',
    'Low-shrub tundra, acidic soils, warmest parts of subzone E':'WBT'
}
Breen_means=Breen_data.groupby('habitat_type').mean().rename(index=ecotype_height_mapping)  
obs_heights['canopy_height_max']=Breen_means['maximum_ canopy_height']
obs_heights['tree_height_mean']=Breen_means['mean_tree_layer_height']
obs_heights['shrub_height_mean']=Breen_means['mean_shrub_layer_height']
obs_heights['tall_shrub_height_mean']=Breen_means['mean_tall_shrub_height']
obs_heights['low_shrub_height_mean']=Breen_means['mean_low_shrub_height']
obs_heights['dwarf_shrub_height_mean']=Breen_means['mean_dwarf_shrub_height']
obs_heights['forb_height_mean']=Breen_means['mean_forb_height']

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

    if not isinstance(moddata['time'].data[0],numpy.datetime64):
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
        for pft in obsdata[landscape_ecotypes[ecotype_num]].index:
            ax.plot([mindate,maxdate],[obsdata[(landscape_ecotypes[ecotype_num],pft)],obsdata[(landscape_ecotypes[ecotype_num],pft)]],c=pft_colors[(pft_names.index(obsdata_PFT_mappings[pft]))],ls='--')

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
    
def plot_PFT_distributions():
    subplot(121)
    names=[]
    bottom=zeros(len(landscape_ecotypes))
    for pftnum in range(len(pft_names_default[:17])):
        pft_pcts=PFT_percents_default.loc[pft_names_default[pftnum]]
        if (pft_pcts==0).all():
            continue
        bar(arange(len(landscape_ecotypes)),pft_pcts,bottom=bottom,facecolor=pft_colors_default[pftnum])
        bottom=bottom+pft_pcts
        names.append(' '.join(pft_names_default[pftnum].split('_')).title() )

    xticks(arange(len(landscape_ecotypes)),landscape_ecotypes,rotation=0)
    title('Default ELM PFTs')
    l=legend(labels=names,loc=(0.0,1.1),fontsize='small')
    l.set_draggable(True)

    subplot(122)
    names=[]
    bottom=zeros(len(landscape_ecotypes))
    for pftnum in range(len(pft_names)):
        pft_pcts=PFT_percents.loc[pft_names[pftnum]]
        if (pft_pcts==0).all():
            continue
        bar(arange(len(landscape_ecotypes)),pft_pcts,bottom=bottom,facecolor=pft_colors[pftnum])
        bottom=bottom+pft_pcts
        names.append(' '.join(pft_names[pftnum].split('_')).title() )

    xticks(arange(len(landscape_ecotypes)),landscape_ecotypes,rotation=0)
    title('Updated ELM PFTs')
    l=legend(labels=names,loc=(0.0,1.1),fontsize='small')
    l.set_draggable(True)

    tight_layout()

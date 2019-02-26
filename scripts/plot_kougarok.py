from pylab import *
import xarray
import numpy
import pandas

params_file='param_files/clm_params_newpfts_c180524_orig.nc'
params_file_default='param_files/clm_params_defaultpfts_c180524_orig.nc'

params_fengming=xarray.open_dataset(params_file,autoclose=True)
params_default=xarray.open_dataset(params_file_default,autoclose=True)
params_new=xarray.open_dataset('clm_params_updated.nc',autoclose=True)

pft_names=[name.strip() for name in params_fengming['pftname'].values.astype(str)]
pft_names_default=[name.strip() for name in params_default['pftname'].values.astype(str)]



# domaindata=xarray.open_dataset('/lustre/or-hydra/cades-ccsi/proj-shared/project_acme/cesminput_ngee/share/domains/domain.clm/domain.lnd.51x63pt_kougarok-NGEE_TransA_navy.nc')
surfdata=xarray.open_dataset('param_files/surfdata_51x63pt_kougarok-NGEE_TransA_simyr1850_c181115-sub12.nc')
surfdata_default=xarray.open_dataset('param_files/surfdata_51x63pt_kougarok-NGEE_TransA_simyr1850_c181115default.nc')

landscape_ecotypes=['NAMC','DSLT','AS','WBT','TTWBT','TT']
ecotype_names={'NAMC':'Non-acidic mountain complex',
               'DSLT':'Dwarf shrub lichen tundra',
               'AS'  :'Alder shrubland',
               'WBT' :'Willow birch tundra',
               'TTWBT':'Tussock tundra/willow birch tundra',
               'TT'  :'Tussock tundra'}

pft_colors=['C%d'%n for n in range(10)] + ['k','purple']

PFT_percents=pandas.DataFrame(data=surfdata.PCT_NAT_PFT.values.squeeze(),index=pft_names,columns=landscape_ecotypes)
PFT_percents_default=pandas.DataFrame(data=surfdata_default.PCT_NAT_PFT.values.squeeze(),index=pft_names_default[:17],columns=landscape_ecotypes)

def read_pftfile(filename,maxyear=None):
    output_PFTs=xarray.open_dataset(filename,autoclose=True)
    if maxyear is not None:
        output_PFTs=output_PFTs.sel(time=array([xx.year for xx in output_PFTs.time.values])<=maxyear)  

    pft_mask=output_PFTs.pfts1d_itype_lunit == 1
    weights=output_PFTs.pfts1d_wtgcell[pft_mask]
    if len(weights) == len(pft_names_default[:17])*len(landscape_ecotypes):
        print('Reading from sim with default PFTs')
        vegdata_PFTs=xarray.Dataset(coords={'time':output_PFTs.time,'PFT':arange(len(pft_names_default[:17])),'ecotype':arange(len(landscape_ecotypes))})
        newshape=(len(landscape_ecotypes),len(pft_names_default[:17]))
    else:
        print('Reading from sim with new PFTs')
        vegdata_PFTs=xarray.Dataset(coords={'time':output_PFTs.time,'PFT':arange(len(pft_names)),'ecotype':arange(len(landscape_ecotypes))})
        newshape=(len(landscape_ecotypes),len(pft_names))

    vegdata_PFTs['weights']=(('ecotype','PFT'),weights.values.reshape(newshape))
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

Koug_meas_biomass=pandas.read_excel('obs_data/NGEEArctic_Q3ELM_KougarokBiomass&NPP_20181112.xlsx',sheet_name='data')\
    .set_index(['Ecotype','ELMgroup'])
Koug_meas_chem=pandas.read_excel('obs_data/NGEEArctic_Q3ELM_KougarokSLA&Chemistry_20181112.xlsx',sheet_name='data')\
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
    
        mindate=minyear
        maxdate=maxyear
        if cumulative:
            dat=dat.cumsum(dim='time')
            dat=dat-dat.isel(time=nonzero(t>minyear)[0][0])
    else:
        t=moddata['time'].data
        mindate=datetime.date(minyear,1,1)
        maxdate=datetime.date(maxyear,1,1)

        if cumulative:
            dat=dat.cumsum(dim='time')
            dat=dat-dat.sel(time=mindate)


    if ax is None:
        ax=gca()

    dat=dat*modfactor
    ax.set_title(longname)
    for pft in moddata.PFT.values:
        pftlabel=pft_names[pft]
        if pftlabel.startswith('arctic_'):
            pftlabel=pftlabel[len('arctic_'):]
        ax.plot(t,dat.sel(ecotype=ecotype_num,PFT=pft),c=pft_colors[pft],label=pftlabel,**kwargs)
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


def save_all_figs(dirname='Figures',format='png',**kwargs):
    for fname in get_figlabels():
        print(fname)
        figure(fname).savefig('{dirname:s}/{fname:s}.{format}'.format(dirname=dirname,format=format,fname=fname),**kwargs)


def pft_params(paramdata,paramnames):
    pftnames=[name.strip() for name in paramdata['pftname'].values.astype(str)]
    if isinstance(paramnames,str):
        paramnames=[paramnames]
    pdict={'pftname':pftnames}
    for name in paramnames:
        pdict[name]=paramdata[name]
    return pandas.DataFrame(pdict).set_index('pftname') 

if __name__=='__main__':
    outputdata_dir='/nfs/data/ccsi/b0u/Kougarok/userpft'

    #vegdata_PFTs_oldparams=read_pftfile(outputdata_dir+'/accelspinup/ELMuserpft_adspinuptest_Kougarok_ICB1850CNPRDCTCBC.h1.nc')
    #vegdata_PFTs_oldparams=read_pftfile(outputdata_dir+'/accelspinup/ELMuserpft_adspinuptest_newparams_Kougarok_ICB1850CNPRDCTCBC.h1.nc',maxyear=150)
    #vegdata_PFTs_newparams=read_pftfile(outputdata_dir+'/accelspinup/ELMuserpft_adspinuptest_newparams_Kougarok_ICB1850CNPRDCTCBC.h1_20190220.nc',maxyear=150)
    vegdata_PFTs_oldparams=read_pftfile(outputdata_dir+'/accelspinup/ELMuserpft_adspinuptest_newparams_Kougarok_ICB1850CNPRDCTCBC.h1_fcur_20190221.nc',maxyear=150)
    #vegdata_PFTs_newparams=read_pftfile(outputdata_dir+'/accelspinup/ELMuserpft_adspinuptest_newparams_Kougarok_ICB1850CNPRDCTCBC.h1_fcur_eg-gram_rhizomefr.nc',maxyear=150)
    #vegdata_PFTs_newparams=read_pftfile(outputdata_dir+'/accelspinup/ELMuserpft_adspinuptest_newparams_Kougarok_ICB1850CNPRDCTCBC.h1_fcur_eg-gram_rhizomefr_leaf-fr-long.nc',maxyear=150)
    vegdata_PFTs_newparams=read_pftfile(outputdata_dir+'/accelspinup/ELMuserpft_adspinuptest_newparams_Kougarok_ICB1850CNPRDCTCBC.h1_2019-02-25.nc')
    #data=xarray.open_dataset(outputdata_dir+'/hist/ELMuserpft_Kougarok_ICB20TRCNPRDCTCBC.clm2.h.nc')
    #data_default=xarray.open_dataset(outputdata_dir+'/hist/ELMuserpft_Kougarok_ICB20TRCNPRDCTCBC_defaultparams.clm2.h.nc')

    minyear=0
    maxyear=50

    meas_leaf_C=(Koug_meas_biomass['LeafBiomass_gperm2']*Koug_meas_chem['LeafC_percent']/100).groupby('Ecotype').sum()
    meas_nonvasc_C=(Koug_meas_biomass['NonvascularBiomass_gperm2']*0.5).groupby('Ecotype').sum()
    # Should this include nonvascular biomass?
    obsdata=meas_leaf_C+meas_nonvasc_C
    # for ecotype in landscape_ecotypes:
    #     obsdata=meas_leaf_C+meas_nonvasc_C
    #     defaultparams_ax.plot([mindate,maxdate],[obsdata[ecotype],obsdata[ecotype]],c='C%d'%(landscape_ecotypes.index(ecotype)),ls=':')
    #     newparams_ax.plot([mindate,maxdate],[obsdata[ecotype],obsdata[ecotype]],c='C%d'%(landscape_ecotypes.index(ecotype)),ls=':')

    meas_root_C=Koug_meas_biomass['FineRootBiomass_gperm2'][:,'mixed']*(Koug_meas_chem['FineRootC_percent']/100).groupby('Ecotype').mean()
    # Should rhizomes be treated as fine or coarse roots?

    meas_stem_C=(Koug_meas_biomass['StemBiomass_gperm2']*Koug_meas_chem['StemC_percent']/100).groupby('Ecotype').sum()
    # Check this. WOODC appears to be sum of LIVESTEMC+DEADSTEMC+LIVECROOTC+DEADCROOTC

    # Try comparing rhizomes to LIVECROOTC
    meas_rhizome_C=(Koug_meas_biomass['RhizomeBiomass_gperm2']*Koug_meas_chem['RhizomeC_percent']/100).groupby('Ecotype').sum()


    figure('PFT distributions');clf()
    subplot(211)
    names=[]
    bottom=zeros(len(landscape_ecotypes))
    for pftnum in range(len(pft_names_default[:17])):
        pft_pcts=PFT_percents_default.loc[pft_names_default[pftnum]]
        if (pft_pcts==0).all():
            continue
        bar(arange(len(landscape_ecotypes)),pft_pcts,bottom=bottom)
        bottom=bottom+pft_pcts
        names.append(pft_names_default[pftnum])

    xticks(arange(len(landscape_ecotypes)),landscape_ecotypes,rotation=0)
    title('Default PFTs')
    legend(labels=names,loc=(1.1,0.0),fontsize='small')

    subplot(212)
    names=[]
    bottom=zeros(len(landscape_ecotypes))
    for pftnum in range(len(pft_names)):
        pft_pcts=PFT_percents.loc[pft_names[pftnum]]
        if (pft_pcts==0).all():
            continue
        bar(arange(len(landscape_ecotypes)),pft_pcts,bottom=bottom,facecolor=pft_colors[pftnum])
        bottom=bottom+pft_pcts
        names.append(pft_names[pftnum])

    xticks(arange(len(landscape_ecotypes)),landscape_ecotypes,rotation=0)
    title('Updated PFTs')
    legend(labels=names,loc=(1.1,0.0),fontsize='small')

    tight_layout()

    meas_leaf_C=(Koug_meas_biomass['LeafBiomass_gperm2']*Koug_meas_chem['LeafC_percent']/100).copy()
    meas_nonvasc_C=(Koug_meas_biomass['NonvascularBiomass_gperm2']*0.5)
    # Should this include nonvascular biomass?
    for ecotype in landscape_ecotypes:
        meas_leaf_C[ecotype,'moss']=meas_nonvasc_C.loc[ecotype,'moss']
        meas_leaf_C[ecotype,'lichen']=meas_nonvasc_C.loc[ecotype,'lichen']
    
    meas_root_C=Koug_meas_biomass['FineRootBiomass_gperm2'][:,'mixed']*(Koug_meas_chem['FineRootC_percent']/100)
    meas_stem_C=(Koug_meas_biomass['StemBiomass_gperm2']*Koug_meas_chem['StemC_percent']/100)
    meas_rhizome_C=(Koug_meas_biomass['RhizomeBiomass_gperm2']*Koug_meas_chem['RhizomeC_percent']/100)
    meas_NPP_C=(Koug_meas_biomass['LeafNPP_gperm2peryr']*Koug_meas_chem['LeafC_percent']/100)+\
        (Koug_meas_biomass['StemNPP_gperm2peryr']*Koug_meas_chem['StemC_percent']/100)#+\

    meas_froot_NPP =(Koug_meas_biomass['FineRootNPP_gperm2peryr']*Koug_meas_chem['FineRootC_percent']/100).groupby('Ecotype').sum()

    plotvars=['leaf','froot','croot','stem','store','npp','cnpp','height']
    nplots=len(plotvars)
    for econum in range(len(landscape_ecotypes)):
        fig=figure(ecotype_names[landscape_ecotypes[econum]],figsize=(15,5))
        fig.clf()

        gs=fig.add_gridspec(ncols=nplots,nrows=2)
        
        subplot_handles={}
        for var in plotvars:
            subplot_handles[var+'_old']=fig.add_subplot(gs[0,plotvars.index(var)])
            subplot_handles[var+'_new']=fig.add_subplot(gs[1,plotvars.index(var)])

        plot_var_PFTs('LEAFC',vegdata_PFTs_oldparams,obsdata=meas_leaf_C,ecotype_num=econum,ax=subplot_handles['leaf_old'])
        plot_var_PFTs('LEAFC',vegdata_PFTs_newparams,obsdata=meas_leaf_C,ecotype_num=econum,ax=subplot_handles['leaf_new'])

        # Should rhizomes be treated as fine or coarse roots?
        plot_var_PFTs('FROOTC',vegdata_PFTs_oldparams,obsdata=meas_root_C,plotsum=True,ecotype_num=econum,ax=subplot_handles['froot_old'])
        plot_var_PFTs('FROOTC',vegdata_PFTs_newparams,obsdata=meas_root_C,plotsum=True,ecotype_num=econum,ax=subplot_handles['froot_new'])

        plot_var_PFTs(['LIVECROOTC','DEADCROOTC'],vegdata_PFTs_oldparams,longname='C Root & rhizome',obsdata=meas_rhizome_C,plotsum=True,ecotype_num=econum,ax=subplot_handles['croot_old'])
        plot_var_PFTs(['LIVECROOTC','DEADCROOTC'],vegdata_PFTs_newparams,longname='C Root & rhizome',obsdata=meas_rhizome_C,plotsum=True,ecotype_num=econum,ax=subplot_handles['croot_new'])

        plot_var_PFTs(['LIVESTEMC','DEADSTEMC'],vegdata_PFTs_oldparams,obsdata=meas_stem_C,longname='Stem C',ecotype_num=econum,ax=subplot_handles['stem_old'])
        plot_var_PFTs(['LIVESTEMC','DEADSTEMC'],vegdata_PFTs_newparams,obsdata=meas_stem_C,longname='Stem C',ecotype_num=econum,ax=subplot_handles['stem_new'])

        plot_var_PFTs('NPP',vegdata_PFTs_oldparams,longname='NPP',modfactor=3600*24,units='gC m$^{-2}$ day $^{-1}$',ecotype_num=econum,ax=subplot_handles['npp_old'])
        plot_var_PFTs('NPP',vegdata_PFTs_newparams,longname='NPP',modfactor=3600*24,units='gC m$^{-2}$ day $^{-1}$',ecotype_num=econum,ax=subplot_handles['npp_new'])

        dat=plot_var_PFTs('AGNPP',vegdata_PFTs_oldparams,longname='AGNPP',modfactor=3600*24,units='gC m$^{-2}$',cumulative=True,obsdata=meas_NPP_C,minyear=maxyear-1,maxyear=maxyear,ecotype_num=econum,ax=subplot_handles['cnpp_old'])
        t=array([tt.year + (tt.dayofyr-1)/365 for tt in dat['time'].data])
        datmax=dat.isel(time=nonzero(t<=maxyear)[0][-1],ecotype=econum).max()
        subplot_handles['cnpp_old'].set_ylim(bottom=-1,top=max(datmax*1.1,meas_NPP_C[landscape_ecotypes[econum]].max()*1.1))

        dat=plot_var_PFTs('AGNPP',vegdata_PFTs_newparams,longname='AGNPP',modfactor=3600*24,units='gC m$^{-2}$',cumulative=True,obsdata=meas_NPP_C,maxyear=maxyear,minyear=maxyear-1,ecotype_num=econum,ax=subplot_handles['cnpp_new'])
        datmax=dat.isel(time=nonzero(t<=maxyear)[0][-1],ecotype=econum).max()
        subplot_handles['cnpp_new'].set_ylim(bottom=-1,top=max(datmax*1.1,meas_NPP_C[landscape_ecotypes[econum]].max()*1.1))

        plot_var_PFTs('HTOP',vegdata_PFTs_oldparams,weight_area=False,ecotype_num=econum,ax=subplot_handles['height_old'])
        plot_var_PFTs('HTOP',vegdata_PFTs_newparams,weight_area=False,ecotype_num=econum,ax=subplot_handles['height_new'])
   
        plot_var_PFTs('STORVEGC',vegdata_PFTs_oldparams,ecotype_num=econum,longname='Stored C',ax=subplot_handles['store_old'])
        plot_var_PFTs('STORVEGC',vegdata_PFTs_newparams,ecotype_num=econum,longname='Stored C',ax=subplot_handles['store_new'])


        subplot_handles['froot_old'].legend(loc=(-1,1.2),ncol=7,fontsize='small')
        
        figtext(0.5,0.95,ecotype_names[landscape_ecotypes[econum]],ha='center',fontsize='large')
        figtext(0.02,0.25,'Updated params',rotation=90,va='center',fontsize='large')
        figtext(0.02,0.75,'Old params',rotation=90,va='center',fontsize='large')

        tight_layout(rect=(0.03,0.0,1.0,0.95))


    # Plot some biomass ratios

    meas_leaf_C=(Koug_meas_biomass['LeafBiomass_gperm2']*Koug_meas_chem['LeafC_percent']/100)
    meas_root_C=(Koug_meas_biomass['FineRootBiomass_gperm2'][:,'mixed']*(Koug_meas_chem['FineRootC_percent']/100))[:,'mixed']
    meas_stem_C=(Koug_meas_biomass['StemBiomass_gperm2']*(Koug_meas_chem['StemC_percent']/100))

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
        #print('{ecotype:s}: {pft:s} leaf frac {leafCfrac:1.2f}, froot_leaf = {froot_leaf:1.2f}'.format(leafCfrac=leafCfrac[ecotype][pft],ecotype=ecotype,pft=pft,froot_leaf=froot_leaf))
        return froot_leaf

    def leaf_stem(ecotype,pft):
        return meas_stem_C[ecotype]/meas_leaf_C[ecotype][pft] 

    leaf=vegdata_PFTs_oldparams['LEAFC_unweighted']
    stemleafratio_oldparams=xarray.Dataset({'stem_leaf_ratio_unweighted':(vegdata_PFTs_oldparams['LIVESTEMC_unweighted']+vegdata_PFTs_oldparams['DEADSTEMC_unweighted'])/leaf.where(leaf>0.01).rolling(time=365*2).max()})
    leaf=vegdata_PFTs_newparams['LEAFC_unweighted']
    stemleafratio_newparams=xarray.Dataset({'stem_leaf_ratio_unweighted':(vegdata_PFTs_newparams['LIVESTEMC_unweighted']+vegdata_PFTs_newparams['DEADSTEMC_unweighted'])/leaf.where(leaf>0.01).rolling(time=365*2).max()})
    frootleafratio_oldparams=xarray.Dataset({'froot_leaf_ratio_unweighted':vegdata_PFTs_oldparams['FROOTC_unweighted']/vegdata_PFTs_oldparams['LEAFC_unweighted']})
    frootleafratio_newparams=xarray.Dataset({'froot_leaf_ratio_unweighted':vegdata_PFTs_newparams['FROOTC_unweighted']/vegdata_PFTs_newparams['LEAFC_unweighted']})

    froot_leaf_obs = meas_root_C/meas_leaf_C.sum(level='Ecotype')

    nplots=2
    for econum in range(len(landscape_ecotypes)):
        fig=figure(ecotype_names[landscape_ecotypes[econum]]+' ratios',figsize=(8,5))
        fig.clf()
        
        frootleaf_old=subplot(2,nplots,1)
        frootleaf_new=subplot(2,nplots,nplots+1)
        stemleaf_old=subplot(2,nplots,2)
        stemleaf_new=subplot(2,nplots,nplots+2)

        plot_var_PFTs('stem_leaf_ratio',stemleafratio_oldparams,obsdata=meas_stem_C/meas_leaf_C,weight_area=False,longname='Stem to leaf ratio',units='gC/gC',ecotype_num=econum,ax=stemleaf_old)
        plot_var_PFTs('stem_leaf_ratio',stemleafratio_newparams,obsdata=meas_stem_C/meas_leaf_C,weight_area=False,longname='Stem to leaf ratio',units='gC/gC',ecotype_num=econum,ax=stemleaf_new)
        stemleaf_new.legend(fontsize='small',ncol=2)
    
        plot_var_PFTs('froot_leaf_ratio',frootleafratio_oldparams,weight_area=False,longname='Froot to leaf ratio',units='gC/gC',ecotype_num=econum,ax=frootleaf_old)
        plot_var_PFTs('froot_leaf_ratio',frootleafratio_newparams,weight_area=False,longname='Froot to leaf ratio',units='gC/gC',ecotype_num=econum,ax=frootleaf_new)
        frootleaf_old.plot([minyear,maxyear],[froot_leaf_obs[landscape_ecotypes[econum]],froot_leaf_obs[landscape_ecotypes[econum]]],'--',c='C0')
        frootleaf_new.plot([minyear,maxyear],[froot_leaf_obs[landscape_ecotypes[econum]],froot_leaf_obs[landscape_ecotypes[econum]]],'--',c='C0')

        #froot_old.legend(loc=(-1,1.2),ncol=7,fontsize='small')
        
        figtext(0.5,0.95,ecotype_names[landscape_ecotypes[econum]]+' ratios',ha='center',fontsize='large')
        figtext(0.02,0.25,'Updated params',rotation=90,va='center',fontsize='large')
        figtext(0.02,0.75,'Old params',rotation=90,va='center',fontsize='large')

        tight_layout(rect=(0.03,0.0,1.0,0.95))
    


    figure('Temperature and root respiration');clf()
    Tsoil10cm=xarray.open_dataset(outputdata_dir+'/accelspinup/ELMuserpft_Kougarok_ICB1850CNPRDCTCBC_clm2_h_20190129.nc',autoclose=True)['TSOI_10CM']
    t2=array([tt.year + (tt.month-.5)/12 for tt in Tsoil10cm.time.data])
    plot(t2,Tsoil10cm.isel(lndgrid=1)-273.15,'b-')
    plot([0,maxyear],[0.0,0.0],'k--')
    ylabel('Soil temperature (C)')

    ax2=twinx()
    plot_var_PFTs('LEAF_MR',vegdata_PFTs_newparams,1,cumulative=False,ax=ax2,modfactor=3600*24,units='gC/m2/day')
    plot_var_PFTs('FROOT_MR',vegdata_PFTs_newparams,1,cumulative=False,ax=ax2,modfactor=3600*24,units='gC/m2/day')
    xlim(10,20)

    tight_layout()


    figure('Temperature and root respiration cumulative');clf()
    Tsoil10cm=xarray.open_dataset(outputdata_dir+'/accelspinup/ELMuserpft_Kougarok_ICB1850CNPRDCTCBC_clm2_h_20190129.nc',autoclose=True)['TSOI_10CM']
    t2=array([tt.year + (tt.month-.5)/12 for tt in Tsoil10cm.time.data])
    plot(t2,Tsoil10cm.isel(lndgrid=1)-273.15,'b-')
    plot([0,maxyear],[0.0,0.0],'k--')
    ylabel('Soil temperature (C)')

    ax2=twinx()
    t=array([tt.year + (tt.dayofyr-1)/365 for tt in vegdata_PFTs_newparams['time'].data])
    ecotype_num=0
    startyear=10
    endyear=25
    for yr in range(startyear,endyear):
        growingseason_start=nonzero((t>=yr)&(vegdata_PFTs_newparams.sel(PFT=3,ecotype=ecotype_num)['GPP_unweighted'].values>0))[0][0]
        growingseason_start_nextyear=nonzero((t>=yr+1)&(vegdata_PFTs_newparams.sel(PFT=3,ecotype=ecotype_num)['GPP_unweighted'].values>0))[0][0]
        xx=arange(growingseason_start,growingseason_start_nextyear)
        plot_var_PFTs('MR',vegdata_PFTs_newparams.isel(time=xx),ecotype_num,cumulative=True,ax=ax2,modfactor=3600*24,units='gC/m2',minyear=yr,ls='--')
        plot_var_PFTs('FROOT_MR',vegdata_PFTs_newparams.isel(time=xx),ecotype_num,cumulative=True,ax=ax2,modfactor=3600*24,units='gC/m2',minyear=yr,ls=':')
        plot_var_PFTs('GPP',vegdata_PFTs_newparams.isel(time=xx),ecotype_num,cumulative=True,ax=ax2,modfactor=3600*24,units='gC/m2')
    xlim(startyear,endyear)

    tight_layout()


    show()

from pylab import *
import xarray
import numpy
import pandas

params_file='param_files/clm_params_newpfts_c180524_orig.nc'
params_file_default='param_files/clm_params_defaultpfts_c180524_orig.nc'

params_data=xarray.open_dataset(params_file)
params_default_data=xarray.open_dataset(params_file_default)

pft_names=[name.strip() for name in params_data['pftname'].values.astype(str)]
pft_names_default=[name.strip() for name in params_default_data['pftname'].values.astype(str)]



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
    output_PFTs=xarray.open_dataset(filename)
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



def plot_var_PFTs(varname,moddata,obsdata=None,minyear=0,maxyear=150,longname=None,units=None,modfactor=1.0,cumulative=False,weight_area=True):
    if isinstance(varname,str):
        dat=moddata[varname+'_unweighted']
        if longname is None:
            longname=dat.long_name
        if units is None:
            units=dat.units
    else:
        # Assuming we are adding multiple data fields together
        dat=moddata[varname[0]+'_unweighted']
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
        t=array([tt.year + (tt.month-.5)/12 for tt in moddata['time'].data])
    
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


    subplot_handles=[]
    for ecotype_num in range(len(landscape_ecotypes)):
        ax=subplot(2,3,ecotype_num+1)
        subplot_handles.append(ax)

        title(ecotype_names[landscape_ecotypes[ecotype_num]])
        for pft in moddata.PFT.values:
            plot(t,dat.sel(ecotype=ecotype_num,PFT=pft)*modfactor,c=pft_colors[pft])
        xlabel('Time')
        ylabel('%s (%s)'%(longname,units))
        ax.set_xlim(left=mindate,right=maxdate)
        xticks(rotation=80)

        if obsdata is not None:
            # Collection was in 2016-2017, but model only goes through 2010. No big deal I guess.
            for pft in obsdata[landscape_ecotypes[ecotype_num]].index:
                plot([mindate,maxdate],[obsdata[(landscape_ecotypes[ecotype_num],pft)],obsdata[(landscape_ecotypes[ecotype_num],pft)]],c=pft_colors[(pft_names.index(obsdata_PFT_mappings[pft]))],ls='--')


    leg_names=[]
    for name in pft_names:
        if name.startswith('arctic_'):
            leg_names.append(name[len('arctic_'):])
        else:
            leg_names.append(name)
    legend(labels=leg_names,fontsize='small',loc=(1.02,0.0))


    tight_layout()
    subplots_adjust(wspace=0.5)
    return subplot_handles


def save_all_figs(dirname='Figures',format='png',**kwargs):
    for fname in get_figlabels():
        print(fname)
        figure(fname).savefig('{dirname:s}/{fname:s}.{format}'.format(dirname=dirname,format=format,fname=fname),**kwargs)


if __name__=='__main__':
    outputdata_dir='/nfs/data/ccsi/b0u/Kougarok/userpft'

    vegdata_PFTs_oldparams=read_pftfile(outputdata_dir+'/accelspinup/ELMuserpft_adspinuptest_Kougarok_ICB1850CNPRDCTCBC.h1.nc')
    vegdata_PFTs_newparams=read_pftfile(outputdata_dir+'/accelspinup/ELMuserpft_adspinuptest_newparams_Kougarok_ICB1850CNPRDCTCBC.h1.nc',maxyear=200)
    #data=xarray.open_dataset(outputdata_dir+'/hist/ELMuserpft_Kougarok_ICB20TRCNPRDCTCBC.clm2.h.nc')
    #data_default=xarray.open_dataset(outputdata_dir+'/hist/ELMuserpft_Kougarok_ICB20TRCNPRDCTCBC_defaultparams.clm2.h.nc')

    minyear=0
    maxyear=150

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

    figure('Leaf biomass (old params)',figsize=(10.2,6.5));clf()
    for ecotype in landscape_ecotypes:
        meas_leaf_C[ecotype,'moss']=meas_nonvasc_C.loc[ecotype,'moss']
        meas_leaf_C[ecotype,'lichen']=meas_nonvasc_C.loc[ecotype,'lichen']
    sp_h=plot_var_PFTs('LEAFC',vegdata_PFTs_oldparams,meas_leaf_C)
    
    figure('Leaf biomass (new params)',figsize=(10.2,6.5));clf()
    for ecotype in landscape_ecotypes:
        meas_leaf_C[ecotype,'moss']=meas_nonvasc_C.loc[ecotype,'moss']
        meas_leaf_C[ecotype,'lichen']=meas_nonvasc_C.loc[ecotype,'lichen']
    sp_h=plot_var_PFTs('LEAFC',vegdata_PFTs_newparams,meas_leaf_C)


    figure('Fine root biomass (old params)',figsize=(10.2,6.5));clf()
    meas_root_C=Koug_meas_biomass['FineRootBiomass_gperm2'][:,'mixed']*(Koug_meas_chem['FineRootC_percent']/100)
    plot_var_PFTs('FROOTC',vegdata_PFTs_oldparams,meas_root_C)
    # Should rhizomes be treated as fine or coarse roots?
    figure('Fine root biomass (new params)',figsize=(10.2,6.5));clf()
    plot_var_PFTs('FROOTC',vegdata_PFTs_newparams,meas_root_C)

    figure('Stem biomass (old params)',figsize=(10.2,6.5));clf()
    meas_stem_C=(Koug_meas_biomass['StemBiomass_gperm2']*Koug_meas_chem['StemC_percent']/100)
    # Check this. WOODC appears to be sum of LIVESTEMC+DEADSTEMC+LIVECROOTC+DEADCROOTC
    plot_var_PFTs(['LIVESTEMC','DEADSTEMC'],vegdata_PFTs_oldparams,meas_stem_C,longname='Stem C')
    figure('Stem biomass (new params)',figsize=(10.2,6.5));clf()
    plot_var_PFTs(['LIVESTEMC','DEADSTEMC'],vegdata_PFTs_newparams,meas_stem_C,longname='Stem C')

    # For now, assuming we can map rhizome biomass to coarse root C. But maybe it should be combined with storage C or something?
    # figure('Coarse root vs rhizome biomass (PFT)',figsize=(10.2,6.5));clf()
    # meas_rhizome_C=(Koug_meas_biomass['RhizomeBiomass_gperm2']*Koug_meas_chem['RhizomeC_percent']/100)
    # plot_var_PFTs(['LIVECROOTC','DEADCROOTC'],meas_rhizome_C,longname='Coarse root C')

    figure('NPP (old params)',figsize=(10.2,6.5));clf()
    # NPP measurements are annual estimates, so probably doesn't make sense to compare with seasonal unless we are plotting cumulative NPP
    # Also need to deal with roots not being species-specific

    # meas_NPP_C=(Koug_meas_biomass['LeafNPP_gperm2peryr']*Koug_meas_chem['LeafC_percent']/100)+\
    #     (Koug_meas_biomass['StemNPP_gperm2peryr']*Koug_meas_chem['StemC_percent']/100)+\
    #     (Koug_meas_biomass['FineRootNPP_gperm2peryr']*Koug_meas_chem['FineRootC_percent']/100)
    plot_var_PFTs('NPP',vegdata_PFTs_oldparams,modfactor=3600*24,units='gC m$^{-2}$ day $^{-1}$')
    figure('NPP (new params)',figsize=(10.2,6.5));clf()
    plot_var_PFTs('NPP',vegdata_PFTs_newparams,modfactor=3600*24,units='gC m$^{-2}$ day $^{-1}$')

    figure('Cumulative AGNPP (old params)',figsize=(10.2,6.5));clf()
    # NPP measurements are annual estimates, so probably doesn't make sense to compare with seasonal unless we are plotting cumulative NPP
    # Also need to deal with roots not being species-specific
    meas_NPP_C=(Koug_meas_biomass['LeafNPP_gperm2peryr']*Koug_meas_chem['LeafC_percent']/100)+\
        (Koug_meas_biomass['StemNPP_gperm2peryr']*Koug_meas_chem['StemC_percent']/100)#+\
    meas_froot_NPP =(Koug_meas_biomass['FineRootNPP_gperm2peryr']*Koug_meas_chem['FineRootC_percent']/100).groupby('Ecotype').sum()

    sp_h=plot_var_PFTs('AGNPP',vegdata_PFTs_oldparams,modfactor=3600*24,units='gC m$^{-2}$',cumulative=True,obsdata=meas_NPP_C,minyear=maxyear-1,maxyear=maxyear)
    for axnum in range(len(sp_h)):
        sp_h[axnum].set_ylim(bottom=-1,top=600)

    figure('Cumulative AGNPP (new params)',figsize=(10.2,6.5));clf()
    sp_h=plot_var_PFTs('AGNPP',vegdata_PFTs_newparams,modfactor=3600*24,units='gC m$^{-2}$',cumulative=True,obsdata=meas_NPP_C,maxyear=maxyear,minyear=maxyear-1)
    for axnum in range(len(sp_h)):
        sp_h[axnum].set_ylim(bottom=-1,top=600)


    #figure('Cumulative BGNPP (PFT)',figsize=(10.2,6.5));clf()
    #meas_rhizomeNPP=(Koug_meas_biomass['RhizomeNPP_gperm2peryr']*Koug_meas_chem['RhizomeC_percent']/100)

    #sp_h=plot_var_PFTs('BGNPP',modfactor=3600*24,units='gC m$^{-2}$',cumulative=True,obsdata=meas_rhizomeNPP,maxyear=2009)

    #total_mod_BGNPP=(vegdata_PFTs['BGNPP_unweighted']*vegdata_PFTs.weights).groupby('ecotype').sum(dim='PFT').cumsum(dim='time')*3600*24
    #total_mod_BGNPP=total_mod_BGNPP-total_mod_BGNPP.sel(time=datetime.date(2008,1,1))
    #for axnum in range(len(sp_h)):
    #    sp_h[axnum].plot(total_mod_BGNPP.time,total_mod_BGNPP.sel(ecotype=axnum) ,'-',c='C0')
    #    obstotal=meas_rhizomeNPP.groupby('Ecotype').sum()[landscape_ecotypes[axnum]]+meas_froot_NPP[landscape_ecotypes[axnum]]
    #    sp_h[axnum].plot([datetime.date(2008,1,1),datetime.date(2009,1,1)],[obstotal,obstotal],'--',c='C0')
    #    sp_h[axnum].set_ylim(bottom=-1,top=600)


    figure('Canopy height (old params)',figsize=(10.2,6.5));clf()
    plot_var_PFTs('HTOP',vegdata_PFTs_oldparams,weight_area=False)
    figure('Canopy height (new params)',figsize=(10.2,6.5));clf()
    plot_var_PFTs('HTOP',vegdata_PFTs_newparams,weight_area=False)
    

    show()


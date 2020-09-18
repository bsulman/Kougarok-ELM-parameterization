from kougarok_plotting import *

Koug_meas_biomass.rename(index={'NAMC':'DLST','TTWBT':'ASV','DSLT':'BEL'},inplace=True)
Koug_meas_chem.rename(index={'NAMC':'DLST','TTWBT':'ASV','DSLT':'BEL'},inplace=True)

meas_leaf_C=(Koug_meas_biomass['Leaf_gperm2']*Koug_meas_chem['LeafC_percent']/100).copy()
meas_nonvasc_C=(Koug_meas_biomass['Nonvascular_gperm2']*0.5)
# Should this include nonvascular biomass?
meas_leaf_C.loc[:,:,'bryophyte'] = meas_nonvasc_C[:,:,'bryophyte']
meas_leaf_C.loc[:,:,'lichen'] = meas_nonvasc_C[:,:,'lichen']

meas_root_C=Koug_meas_biomass['FineRoot_gperm2'][:,:,'mixed']*(Koug_meas_chem['FineRootC_percent'][:,:,'mixed']/100)
meas_stem_C=(Koug_meas_biomass['Stem_gperm2']*Koug_meas_chem['StemC_percent']/100)
meas_rhizome_C=(Koug_meas_biomass['Rhizome_gperm2']*Koug_meas_chem['RhizomeC_percent']/100)
meas_NPP_C=(Koug_meas_biomass['LeafNPP_gperm2peryear']*Koug_meas_chem['LeafC_percent']/100)+\
    (Koug_meas_biomass['StemNPP_gperm2peryear']*Koug_meas_chem['StemC_percent']/100)#+\

meas_froot_NPP =(Koug_meas_biomass['FineRootNPP_gperm2peryear']*Koug_meas_chem['FineRootC_percent']/100)[:,:,'mixed']
meas_leaf_NPP =(Koug_meas_biomass['LeafNPP_gperm2peryear']*Koug_meas_chem['LeafC_percent']/100)
meas_stem_NPP =(Koug_meas_biomass['StemNPP_gperm2peryear']*Koug_meas_chem['StemC_percent']/100)
meas_rhizome_NPP =(Koug_meas_biomass['RhizomeNPP_gperm2peryear']*Koug_meas_chem['RhizomeC_percent']/100)

# Add together willow and birch tall shrubs
meas_leaf_C.loc[:,:,'potential tall shrub deciduous birch'] = meas_leaf_C.loc[:,:,'potential tall shrub deciduous birch']+meas_leaf_C.loc[:,:,'potential tall shrub deciduous willow'] 
meas_leaf_C.rename(index={'potential tall shrub deciduous birch':'potential tall shrub deciduous non-alder'},inplace=True)
meas_leaf_C.drop('potential tall shrub deciduous willow',level='ELM_PFT',inplace=True)
meas_stem_C.loc[:,:,'potential tall shrub deciduous birch'] = meas_stem_C.loc[:,:,'potential tall shrub deciduous birch']+meas_stem_C.loc[:,:,'potential tall shrub deciduous willow'] 
meas_stem_C.rename(index={'potential tall shrub deciduous birch':'potential tall shrub deciduous non-alder'},inplace=True)
meas_stem_C.drop('potential tall shrub deciduous willow',level='ELM_PFT',inplace=True)
meas_rhizome_C.loc[:,:,'potential tall shrub deciduous birch'] = meas_rhizome_C.loc[:,:,'potential tall shrub deciduous birch']+meas_rhizome_C.loc[:,:,'potential tall shrub deciduous willow'] 
meas_rhizome_C.rename(index={'potential tall shrub deciduous birch':'potential tall shrub deciduous non-alder'},inplace=True)
meas_rhizome_C.drop('potential tall shrub deciduous willow',level='ELM_PFT',inplace=True)
meas_NPP_C.loc[:,:,'potential tall shrub deciduous birch'] = meas_NPP_C.loc[:,:,'potential tall shrub deciduous birch']+meas_NPP_C.loc[:,:,'potential tall shrub deciduous willow'] 
meas_NPP_C.rename(index={'potential tall shrub deciduous birch':'potential tall shrub deciduous non-alder'},inplace=True)
meas_NPP_C.drop('potential tall shrub deciduous willow',level='ELM_PFT',inplace=True)

Amy_Verity_mapping={'deciduous dwarf shrub':'dwarf shrub deciduous',
             'evergreen shrub':'dwarf shrub evergreen',
             'deciduous low shrub':'low shrub deciduous',
             'deciduous low to tall willow and birch shrub':'potential tall shrub deciduous non-alder',
             'deciduous low to tall alder shrub':'potential tall shrub deciduous alder'}
pftareas=pandas.read_excel('../../Amy_data/ngee_arctic_pft_tallies_kougarok_20191010.xlsx',sheet_name='KG_species_table',skiprows=[1,2],header=0).replace(Amy_Verity_mapping) 

pftareas['ELM_default_PFT']=[obsdata_E3SMPFT_mappings[pft] for pft in pftareas['NGEE_PFT']] 

def get_taxon_areas(data,pft,community,pftcol='NGEE_PFT'):

    com_cols=['KG_%s%d_VgComp'%(community,num) for num in range(1,6)]
    if pft not in data[pftcol].to_list():
        raise ValueError('PFT %s not in dataset'%pft)
    return data[data[pftcol]==Amy_Verity_mapping.get(pft,pft)][['Origina_ taxon_determination']+com_cols]

# # Reclassify Betula nana in WBT and ASV as low shrub (which is Amy's classification)
# # Probably need to revisit other communities to make sure Betula nana isn't causing problems
# Update (May 4, 2020): Verity has updated the data to correct the classification issues
# for comm in ['WBT','ASV','TT']:
#     for plotid in unique(meas_leaf_C.loc[comm].index.get_level_values('PlotID')): # Indexing ['WBT',:,pft] didn't work right
#         meas_leaf_C.loc[comm,plotid,'low shrub deciduous'] = meas_leaf_C.loc[comm,plotid,'low shrub deciduous'] + meas_leaf_C.loc[comm,plotid,'dwarf shrub deciduous'] 
#         meas_leaf_C.loc[comm,plotid,'dwarf shrub deciduous']=0.0
#         meas_stem_C.loc[comm,plotid,'low shrub deciduous'] = meas_stem_C.loc[comm,plotid,'low shrub deciduous'] + meas_stem_C.loc[comm,plotid,'dwarf shrub deciduous'] 
#         meas_stem_C.loc[comm,plotid,'dwarf shrub deciduous']=0.0
#         meas_rhizome_C.loc[comm,plotid,'low shrub deciduous'] = meas_rhizome_C.loc[comm,plotid,'low shrub deciduous'] + meas_rhizome_C.loc[comm,plotid,'dwarf shrub deciduous'] 
#         meas_rhizome_C.loc[comm,plotid,'dwarf shrub deciduous']=0.0
#         meas_NPP_C.loc[comm,plotid,'low shrub deciduous'] = meas_NPP_C.loc[comm,plotid,'low shrub deciduous'] + meas_NPP_C.loc[comm,plotid,'dwarf shrub deciduous'] 
#         meas_NPP_C.loc[comm,plotid,'dwarf shrub deciduous']=0.0
# 



def obs_to_defaultPFTs(obsdata,mapping=obsdata_defaultPFT_mappings):
    # Best way I could figure out how to do this, sorry. Making a multilevel index with the same community and plot levels but new PFTs
    new_index=pandas.MultiIndex.from_tuples([(x[0][0],x[0][1],x[1]) for x in pandas.MultiIndex.from_product((obsdata.index.droplevel('ELM_PFT').drop_duplicates().values,unique(list(mapping.values())))).values],names=('Ecotype','PlotID','ELM_PFT')) 
    if isinstance(obsdata,pandas.Series):
        out=pandas.Series(index=new_index,data=0)
        for pft in unique(obsdata.index.get_level_values('ELM_PFT')):
            out.loc[:,:,mapping[pft]] = out.loc[:,:,mapping[pft]].add( obsdata.loc[:,:,pft].fillna(0.0), fill_value=0.0)

    else:
        raise TypeError('Must be a Series')

    return out
    
def new_to_defaultPFTs(data,container_template,fields=('LEAFC','FROOTC','LIVECROOTC','DEADCROOTC','STORVEGC','LIVESTEMC','DEADSTEMC',)):
    dataset_out=container_template.coords.to_dataset()
    shrub=container_template.PFTnames.data.tolist().index('broadleaf_deciduous_boreal_shrub')  
    shrub_new=nonzero(data.PFTnames.str.contains('shrub').data)[0]  
    grass=container_template.PFTnames.data.tolist().index('c3_arctic_grass')
    grass_new=nonzero(data.PFTnames.str.contains('graminoid').data|data.PFTnames.str.contains('forb').data)[0]
    nonveg=0
    nonveg_new=array([0,1,2]) # Not vegetated and two nonvascular plant types. But this is inconsistent with model level 3, which treats nonvascular as grasses
    
    for field in fields:
        # These should be weighted averages
        container=xarray.DataArray(coords=container_template[field+'_unweighted'].coords)
        v=container.isel(PFT=shrub)
        v[:]=(data[field+'_unweighted'].isel(PFT=shrub_new)*data['weights'].isel(PFT=shrub_new)).sum(dim='PFT')/data['weights'].isel(PFT=shrub_new).sum(dim='PFT')
        # for pft in shrub_new:
        #     container[shrub,:,:] = container[shrub,:,:]+data[field+'_unweighted'][:,:,pft]
        v=container.isel(PFT=grass)
        v[:]=(data[field+'_unweighted'].isel(PFT=grass_new)*data['weights'].isel(PFT=grass_new)).sum(dim='PFT')/data['weights'].isel(PFT=grass_new).sum(dim='PFT')
        # for pft in grass_new:
        #     container[grass,:,:] = container[grass,:,:]+data[field+'_unweighted'][:,:,pft]
        v=container.isel(PFT=nonveg)
        v[:]=(data[field+'_unweighted'].isel(PFT=nonveg_new)*data['weights'].isel(PFT=nonveg_new)).sum(dim='PFT')/data['weights'].isel(PFT=nonveg_new).sum(dim='PFT')
        # for pft in nonveg_new:
        #     container[nonveg,:,:] = container[nonveg,:,:]+data[field+'_unweighted'][:,:,pft]
            
        dataset_out[field+'_unweighted']=container
        
    return dataset_out
    
def plot_mod_bar_stack(x,dat,econum,do_legend=False,op='max',**kwargs):
    handles=[]
    bottom=0.0
    for pftnum in dat.PFT:
        if dat.PFTnames.sel(PFT=pftnum)=='arctic_dry_graminoid':
            val=getattr(dat.sel(PFT=pftnum,ecotype=econum),op)(dim='time')+getattr(dat.sel(PFT=pftnum+1,ecotype=econum).fillna(0.0),op)(dim='time')
        elif dat.PFTnames.sel(PFT=pftnum)=='arctic_wet_graminoid':
            continue
        else:
            val=getattr(dat.sel(PFT=pftnum,ecotype=econum),op)(dim='time')
        if ~isnan(val):
            name=prettify_pft_name(dat.PFTnames.sel(PFT=pftnum).item())
            handles.append(bar(x,val,bottom=bottom,facecolor=dat.PFTcolors.sel(PFT=pftnum).item(),label=name,**kwargs))
            bottom=bottom+val
    return handles


def plot_obs_bar_stack(x,obsdata,ecotype_num,bottom=0.0,pfts=None,**kwargs):
    names=[]
    if pfts is None:
        pfts=unique(obsdata[landscape_ecotypes[ecotype_num]].index.get_level_values('ELM_PFT'))
    for num,pft in enumerate(pfts):
        if 'PlotID' in obsdata.index.names:
            val=obsdata.loc[landscape_ecotypes[ecotype_num],:,pft].mean()
            valstd=obsdata.loc[landscape_ecotypes[ecotype_num],:,pft].std()
        else:
            val=obsdata[(landscape_ecotypes[ecotype_num],pft)]
            valstd=nan

        if ~isnan(val):
            bar(x,val,bottom=bottom,facecolor=pft_colors[(pft_names.index(obsdata_PFT_mappings[pft]))],**kwargs)
            errorbar(x+num/len(pfts)*0.25*kwargs.get('width',1.0),bottom+val,yerr=valstd,fmt='none',color=pft_colors[(pft_names.index(obsdata_PFT_mappings[pft]))])
            bottom=bottom+val


def plot_pft_biomass_bars(data_to_plot,axs=None,per_area=False,include_froot=True,include_storage=True,use_pooled_default_PFTs=True,pfts_inuse=None,leg_axis=0):
    if pfts_inuse is None:
        pfts_inuse=data_to_plot['PFT'][(data_to_plot.weights>0).any(dim='ecotype')].data.tolist()
    w=0.8/6/2
    if 'c3_arctic_grass' in data_to_plot.PFTnames and use_pooled_default_PFTs:
        if axs is None:
            ncols=len(pfts_inuse)
            nrows=1
            axs=gcf().subplots(ncols=ncols,nrows=nrows,squeeze=False).ravel()
        for n in range(len(pfts_inuse[1:])):
            axs[n+1].set_title(prettify_pft_name(pft_names_default[pfts_inuse[n+1]]))
            axs[n+1].set_ylabel('Biomass (gC m$^{-2}$)')
            axs[n+1].set_xticks(arange(6)/6)
            axs[n+1].set_xticklabels(landscape_ecotypes)
        sca(axs[0])
        title('Nonvascular')
        ylabel('Biomass (gC m$^{-2}$)')
        xticks(arange(6)/6,landscape_ecotypes)
    else:
        if axs is None:
            nrows=2
            ncols=5
            axs=gcf().subplots(ncols=ncols,nrows=nrows).ravel()
        for n in range(len(pfts_inuse[1:])):
            sca(axs[n])
            title(prettify_pft_name(pft_names[pfts_inuse[n+1]]))
            ylabel('Biomass (gC m$^{-2}$)')
            xticks(arange(6)/6,landscape_ecotypes)
    for pftnum in pfts_inuse:
        if 'c3_arctic_grass' in data_to_plot.PFTnames:
            if use_pooled_default_PFTs:
                n=pfts_inuse.index(pftnum)+1
            else:
                n=pft_names.index(default_new_mappings[data_to_plot.PFTnames.values[pftnum]])

        else:
            n=pfts_inuse.index(pftnum)

        sca(axs[n-1])
        bottom=zeros(6)
        bottom_obs=zeros(6)
        handles=[]
        
        if include_froot:
            vals=get_var_PFTs('FROOTC',data_to_plot,weight_area=not per_area).max(dim='time').sel(PFT=pftnum)
            h=bar(arange(6)/6,vals,bottom=bottom,width=w,facecolor='brown',label='Fine root',hatch='//')
            bottom=bottom+vals.fillna(0).data
            handles.append(h)
        
        vals=get_var_PFTs(['LIVECROOTC','DEADCROOTC'],data_to_plot,weight_area=not per_area).mean(dim='time').sel(PFT=pftnum)
        h=bar(arange(6)/6,vals,bottom=bottom,width=w,facecolor='orange',label='Coarse root',hatch='//')
        bottom=bottom+vals.fillna(0).data
        handles.append(h)
        
        if include_storage:
            vals=get_var_PFTs('STORVEGC',data_to_plot,weight_area=not per_area).min(dim='time').sel(PFT=pftnum)
            h=bar(arange(6)/6,vals,bottom=bottom,width=w,facecolor='yellow',label='Storage',hatch='//')
            bottom=bottom+vals.fillna(0).data
            handles.append(h)
            # bar(arange(6)/6+w,meas_rhizome_C)
        
        vals=get_var_PFTs(['LIVESTEMC','DEADSTEMC'],data_to_plot,weight_area=not per_area).mean(dim='time').sel(PFT=pftnum)
        h=bar(arange(6)/6,vals,bottom=bottom,width=w,facecolor='blue',label='Stem',hatch='//')
        bottom=bottom+vals.fillna(0).data
        handles.append(h)
        
        vals=get_var_PFTs('LEAFC',data_to_plot,weight_area=not per_area).max(dim='time').sel(PFT=pftnum)
        h=bar(arange(6)/6,vals,bottom=bottom,width=w,facecolor='green',label='Leaf',hatch='//')
        bottom=bottom+vals.fillna(0).data
        handles.append(h)  

        
    
    if 'c3_arctic_grass' in data_to_plot.PFTnames and use_pooled_default_PFTs:
        if 9 not in pfts_inuse:
            mapping=obsdata_E3SMPFT_mappings
        else:
            mapping=obsdata_defaultPFT_mappings
        obsdata={'leaf':obs_to_defaultPFTs(meas_leaf_C,mapping),
                 'rhizome':obs_to_defaultPFTs(meas_rhizome_C,mapping),
                 'stem':obs_to_defaultPFTs(meas_stem_C,mapping),
                 'nonvasc':obs_to_defaultPFTs(meas_nonvasc_C,mapping)}
        area_taxon_col='ELM_default_PFT'
    else:
        obsdata={'leaf':meas_leaf_C,
                 'rhizome':meas_rhizome_C,
                 'stem':meas_stem_C,
                 'nonvasc':meas_nonvasc_C}                
        area_taxon_col='NGEE_PFT'
    for pftname in unique(obsdata['leaf'].index.get_level_values('ELM_PFT')) :
        if pftname in ['other','mixed','not_vegetated']:
             continue
        if 'c3_arctic_grass' in data_to_plot.PFTnames and use_pooled_default_PFTs:
            if pftname=='nonvascular':
                ax=axs[0]
            else:
                ax=axs[pfts_inuse.index(data_to_plot.PFTnames.values.tolist().index(pftname))]
        else:
            pftnum=pft_names.index(obsdata_PFT_mappings[pftname])
            if pftnum not in pfts_inuse:
                continue
            ax=axs[pfts_inuse.index(pftnum)-1]
        for econum in range(len(landscape_ecotypes)):
            # if per_area:
            #     if 'c3_arctic_grass' in data_to_plot.PFTnames and use_pooled_default_PFTs:
            #         if pftname=='nonvascular':
            #             areafrac=(PFT_percents[landscape_ecotypes[econum]]['arctic_lichen']+PFT_percents[landscape_ecotypes[econum]]['arctic_bryophyte'])/100
            #         else:
            #             areafrac=PFT_percents_default[landscape_ecotypes[econum]][pftname]/100
            #     else:
            #         areafrac=PFT_percents[landscape_ecotypes[econum]][obsdata_PFT_mappings[pftname]]/100
            #     if areafrac==0.0:
            #         invareafrac=0.0
            #     else:
            #         invareafrac=1.0/areafrac
            # else:
            #     invareafrac=1.0
            plotIDs=unique(obsdata['leaf'][landscape_ecotypes[econum]].index.get_level_values('PlotID'))
            for plotnum,plotID in enumerate(plotIDs):
                bottom=0.0
                if per_area:
                    areafrac=get_taxon_areas(pftareas,pftname,landscape_ecotypes[econum],area_taxon_col)['KG_%s%s_VgComp'%(landscape_ecotypes[econum],plotID.split('-')[-1])].sum()/100
                else:
                    areafrac=1.0
                if areafrac==0.0:
                    invareafrac=0.0
                else:
                    invareafrac=1.0/areafrac
                if per_area and pftname not in ['moss','bryophyte','lichen','nonvascular'] and pftname in obsdata['leaf'][landscape_ecotypes[econum],plotID] and obsdata['leaf'][landscape_ecotypes[econum],plotID,pftname]>0:
                    # print(pftname,plotID,invareafrac,obsdata['leaf'][landscape_ecotypes[econum],plotID,pftname])
                    val=meas_root_C[landscape_ecotypes[econum],plotID] # Does this need to be adjusted by fractional areas at all?
                    h_obsroot=ax.bar(econum/6+w+plotnum*w/len(plotIDs),val,bottom=bottom,facecolor='brown',width=w/len(plotIDs),label='Fine roots')
                    bottom=bottom+val
                
                val=obsdata['rhizome'][landscape_ecotypes[econum],plotID].fillna(0.0).get(pftname,0.0)*invareafrac
                ax.bar(econum/6+w+plotnum*w/len(plotIDs),val,bottom=bottom,facecolor='orange',width=w/len(plotIDs))
                bottom=bottom+val
                val=obsdata['stem'][landscape_ecotypes[econum],plotID].fillna(0.0).get(pftname,0.0)*invareafrac
                ax.bar(econum/6+w+plotnum*w/len(plotIDs),val,bottom=bottom,facecolor='blue',width=w/len(plotIDs))
                bottom=bottom+val
                val=obsdata['leaf'][landscape_ecotypes[econum],plotID].fillna(0.0).get(pftname,0.0)*invareafrac
                ax.bar(econum/6+w+plotnum*w/len(plotIDs),val,bottom=bottom,facecolor='green',width=w/len(plotIDs))
                bottom=bottom+val
                val=obsdata['nonvasc'][landscape_ecotypes[econum],plotID].fillna(0.0).get(pftname,0.0)*invareafrac
                h_nonvasc=ax.bar(econum/6+w+plotnum*w/len(plotIDs),val,bottom=bottom,facecolor='cyan',width=w/len(plotIDs),label='Non-vascular')
                bottom=bottom+val
                


    handles.append(h_nonvasc)
    if leg_axis is not None:
        l=axs[leg_axis].legend(fontsize='medium',ncol=1,handles=handles)
        l.set_draggable(True)



if __name__=='__main__':
    
    import warnings

    # Don't need to print this warning a billion times
    warnings.filterwarnings(action='ignore',message='All-NaN slice encountered')  
    warnings.filterwarnings(action='ignore',message='Mean of empty slice')  

    from argparse import ArgumentParser
    
    parser = ArgumentParser()
    parser.add_argument('filename',help='netCDF file to plot')
    parser.add_argument('--all_timeseries',action='store_true')
    parser.add_argument('--cumul_rootresp',action='store_true')
    parser.add_argument('--SOMstocks',action='store_true')
    parser.add_argument('--Nfix',action='store_true')
    parser.add_argument('--PFTdist',action='store_true')
    parser.add_argument('--biomass',action='store_true')
    parser.add_argument('--height',action='store_true')
    parser.add_argument('--AG_BG',action='store_true')
    parser.add_argument('--nonvasc',action='store_true')
    parser.add_argument('--vcmax',action='store_true')
    parser.add_argument('--sitemap',action='store_true')
    parser.add_argument('--rootprofile',action='store_true')
    parser.add_argument('--all',action='store_true')
    parser.add_argument('--save_figs',help='Directory to save figures')

    options = parser.parse_args()

    filename=options.filename
        
    import os
    dataname=os.path.basename(filename)
    vegdata_PFTs=read_pftfile(filename,maxyear=None)
    
    columndata=xarray.open_dataset(filename.replace('h1','h0').replace('h2','h0'))
    t_col=array([tt.year + (tt.dayofyr-1)/365 for tt in columndata['time'].data])

    minyear=int(floor(t_col.min()))
    maxyear=int(ceil(t_col.max()))
    


    plotvars=['leaf','froot','croot','store','stem','cnpp','height','downreg']
    nplots=len(plotvars)
    
    # if 'c3_arctic_grass' in vegdata_PFTs.PFTnames:
    #     pfts_inuse=[]
    #     for pftnum in range(len(vegdata_PFTs.PFTnames)):
    #         pft_pcts=PFT_percents_default.loc[vegdata_PFTs.PFTnames[pftnum].values]
    #         if not (pft_pcts==0).all():
    #             pfts_inuse.append(pftnum)
    # else:
    #     pfts_inuse=vegdata_PFTs.PFT.values[:-1]
    pfts_inuse=vegdata_PFTs['PFT'][(vegdata_PFTs.weights>0).any(dim='ecotype')].data.tolist()
        
        
    year=xarray.DataArray(array([tt.year for tt in vegdata_PFTs.time.data]),dims='time',coords={'time':vegdata_PFTs.time},name='year')
    annualmax=vegdata_PFTs.groupby(year,squeeze=False).max(dim='time',keep_attrs=True)
    annualmax['year']=vegdata_PFTs['time'].groupby(year).first()
    annualmax=annualmax.rename({'year':'time'})
    annualmean=vegdata_PFTs.groupby(year,squeeze=False).mean(dim='time',keep_attrs=True)
    annualmean['year']=vegdata_PFTs['time'].groupby(year).first()
    annualmean=annualmean.rename({'year':'time'})
        
    if options.all_timeseries or options.all:
        
        fig=figure('Ecotype time series (%s)'%dataname,figsize=(15,10))
        fig.clf()
        gs=fig.add_gridspec(ncols=nplots,nrows=len(vegdata_PFTs.ecotype))
        subplot_handles={}
        
        # demand_ratio = xarray.Dataset({'Nratio_unweighted':vegdata_PFTs['SMINN_TO_NPOOL_unweighted']/vegdata_PFTs['PLANT_NDEMAND_unweighted'],
        #                                 'Pratio_unweighted':vegdata_PFTs['SMINP_TO_PPOOL_unweighted']/vegdata_PFTs['PLANT_PDEMAND_unweighted']})
        # demand_ratio_annual = xarray.Dataset({'Nratio_unweighted':annualmean['SMINN_TO_NPOOL_unweighted']/annualmean['PLANT_NDEMAND_unweighted'],
        #                                 'Pratio_unweighted':annualmean['SMINP_TO_PPOOL_unweighted']/annualmean['PLANT_PDEMAND_unweighted']})

        for econum in range(len(vegdata_PFTs.ecotype)):
            ecotype=landscape_ecotypes[econum]
            subplot_handles[ecotype]={}
            for var in plotvars:
                subplot_handles[ecotype][var]=fig.add_subplot(gs[econum,plotvars.index(var)])

            plot_var_PFTs('LEAFC',annualmax,obsdata=meas_leaf_C,ecotype_num=econum,ax=subplot_handles[ecotype]['leaf'],minyear=minyear,maxyear=maxyear)
            subplot_handles[ecotype]['leaf'].text(-1.0,0.5,ecotype,rotation=90,va='center',ha='center',transform=subplot_handles[ecotype]['leaf'].transAxes)
            # Should rhizomes be treated as fine or coarse roots?
            plot_var_PFTs('FROOTC',annualmax,obsdata=meas_root_C,plotsum=True,ecotype_num=econum,ax=subplot_handles[ecotype]['froot'],minyear=minyear,maxyear=maxyear)

            plot_var_PFTs(['LIVECROOTC','DEADCROOTC'],vegdata_PFTs,longname='C Root',plotsum=True,ecotype_num=econum,ax=subplot_handles[ecotype]['croot'],minyear=minyear,maxyear=maxyear)

            plot_var_PFTs(['LIVESTEMC','DEADSTEMC'],vegdata_PFTs,obsdata=meas_stem_C,longname='Stem C',ecotype_num=econum,ax=subplot_handles[ecotype]['stem'],minyear=minyear,maxyear=maxyear)

            # plot_var_PFTs('NPP',vegdata_PFTs,longname='NPP',modfactor=3600*24,units='gC m$^{-2}$ day $^{-1}$',ecotype_num=econum,ax=subplot_handles[ecotype]['npp'],maxyear=maxyear)

            dat=plot_var_PFTs('AGNPP',vegdata_PFTs,longname='AGNPP',modfactor=3600*24*365,units='gC m$^{-2}$',cumulative=True,obsdata=meas_NPP_C,maxyear=maxyear,minyear=maxyear-1,ecotype_num=econum,ax=subplot_handles[ecotype]['cnpp'])
            t=array([tt.year + (tt.dayofyr-1)/365 for tt in dat['time'].data])
            datmax=dat.isel(time=nonzero(t<=maxyear)[0][-1],ecotype=econum).max()
            subplot_handles[ecotype]['cnpp'].set_ylim(bottom=-1,top=max(datmax*1.1,meas_NPP_C[ecotype].max()*1.1))

            plot_var_PFTs('HTOP',vegdata_PFTs,weight_area=False,ecotype_num=econum,ax=subplot_handles[ecotype]['height'],minyear=minyear,maxyear=maxyear)

            plot_var_PFTs('STORVEGC',vegdata_PFTs,ecotype_num=econum,longname='Stored C',obsdata=meas_rhizome_C,ax=subplot_handles[ecotype]['store'],minyear=minyear,maxyear=maxyear)

            
            # plot_var_PFTs('DOWNREG',annualmax,weight_area=False,ecotype_num=econum,ax=subplot_handles[ecotype]['downreg'],maxyear=maxyear,longname='Nut. lim',units='frac')
            subplot_handles[ecotype]['downreg'].plot(((columndata['NPP']*columndata['FPG']).groupby(year).sum(dim='time')/columndata['NPP'].groupby(year).sum(dim='time')).isel(lndgrid=econum))
            subplot_handles[ecotype]['downreg'].plot(((columndata['NPP']*columndata['FPG_P']).groupby(year).sum(dim='time')/columndata['NPP'].groupby(year).sum(dim='time')).isel(lndgrid=econum))
            title('Nutrient lim.')

            for ax in subplot_handles[ecotype].values():
                ax.title.set_fontsize('small')
                ax.xaxis.label.set_fontsize('small')
                ax.yaxis.label.set_fontsize('small')
                for label in ax.get_xticklabels():
                    label.set_fontsize('small')
                for label in ax.get_yticklabels():
                    label.set_fontsize('small')

        

        tight_layout(rect=(0.03,0.0,1.0,0.95))
        subplots_adjust(hspace=0.8)
        subplot_handles[landscape_ecotypes[0]]['froot'].legend(loc=(-1,1.2),ncol=7,fontsize='small',handles=list(array(subplot_handles[landscape_ecotypes[0]]['froot'].lines)[pfts_inuse]))

    if options.cumul_rootresp or options.all:
        Tsoil10cm=columndata['TSOI_10CM']

        figure('Temperature and root respiration cumulative (%s)'%dataname);clf()
        
        t2=array([tt.year + (tt.month-.5)/12 for tt in Tsoil10cm.time.data])
        plot(t2,Tsoil10cm.isel(lndgrid=min(1,len(Tsoil10cm.lndgrid)-1))-273.15,'b-')
        plot([0,maxyear],[0.0,0.0],'k--')
        ylabel('Soil temperature (C)')

        ax2=twinx()
        t=array([tt.year + (tt.dayofyr-1)/365 for tt in vegdata_PFTs['time'].data])
        ecotype_num=min(1,len(vegdata_PFTs.ecotype)-1)
        startyear=minyear+10
        endyear=minyear+25
        for yr in range(startyear,endyear):
            growingseason_start=nonzero((t>=yr)&(vegdata_PFTs.sel(ecotype=ecotype_num)['GPP_unweighted'].values>0).any(axis=1))[0][0]
            growingseason_start_nextyear=nonzero((t>=yr+1)&(vegdata_PFTs.sel(ecotype=ecotype_num)['GPP_unweighted'].values>0).any(axis=1))[0][0]
            xx=arange(growingseason_start,growingseason_start_nextyear)
            plot_var_PFTs('MR',vegdata_PFTs.isel(time=xx),ecotype_num,cumulative=True,ax=ax2,modfactor=3600*24,units='gC/m2',minyear=yr,ls='--')
            plot_var_PFTs('FROOT_MR',vegdata_PFTs.isel(time=xx),ecotype_num,cumulative=True,ax=ax2,modfactor=3600*24,units='gC/m2',minyear=yr,ls=':')
            plot_var_PFTs('LIVECROOT_MR',vegdata_PFTs.isel(time=xx),ecotype_num,cumulative=True,ax=ax2,modfactor=3600*24,units='gC/m2',minyear=yr,ls=':',lw=0.5)
            plot_var_PFTs('LEAF_MR',vegdata_PFTs.isel(time=xx),ecotype_num,cumulative=True,ax=ax2,modfactor=3600*24,units='gC/m2',minyear=yr,ls='-.')
            plot_var_PFTs('GPP',vegdata_PFTs.isel(time=xx),ecotype_num,cumulative=True,ax=ax2,modfactor=3600*24,units='gC/m2')
        xlim(startyear,endyear)

        tight_layout()
    
    if options.SOMstocks or options.all:
        
        barplots=True
        figure('SOM stocks (%s)'%dataname,figsize=(10,5));clf()
        if not barplots:
            subplot(131)
            plot(t_col,columndata.TOTSOMC/1000)
            title('Soil C')
            xlabel('Time (years)')
            ylabel('{name:s} ({units:s})'.format(name=columndata.TOTSOMC.long_name,units='kgC m$^{-2}$'))
            subplot(132)
            plot(t_col,columndata.TOTSOMN/1000)
            title('Soil N')
            xlabel('Time (years)')
            ylabel('{name:s} ({units:s})'.format(name=columndata.TOTSOMN.long_name,units='kgN m$^{-2}$'))
            subplot(133)
            plot(t_col,columndata.TOTSOMP)
            title('Soil P')
            xlabel('Time (years)')
            ylabel('{name:s} ({units:s})'.format(name=columndata.TOTSOMP.long_name,units=columndata.TOTSOMP.units))
            

        else:
            x=arange(6)
            subplot(131)
            bar(x,columndata.TOTSOMC.sel(time=t_col>(t_col.max()-100)).mean(dim='time')/1000,hatch='//')
            title('Soil C')
            xticks(x,landscape_ecotypes)
            ylabel('{name:s} ({units:s})'.format(name=columndata.TOTSOMC.long_name,units='kgC m$^{-2}$'))
            subplot(132)
            bar(x,columndata.TOTSOMN.sel(time=t_col>(t_col.max()-100)).mean(dim='time')/1000,hatch='//')
            title('Soil N')
            xticks(x,landscape_ecotypes)
            ylabel('{name:s} ({units:s})'.format(name=columndata.TOTSOMN.long_name,units='kgN m$^{-2}$'))
            subplot(133)
            bar(x,columndata.TOTSOMP.sel(time=t_col>(t_col.max()-100)).mean(dim='time'),hatch='//')
            title('Soil P')
            xticks(x,landscape_ecotypes)
            ylabel('{name:s} ({units:s})'.format(name=columndata.TOTSOMP.long_name,units=columndata.TOTSOMP.units))
            

        tight_layout()

    if options.Nfix or options.all:
        
        figure('N fixation(%s)'%dataname,figsize=(6.4,6.7));clf()
        subplot(211)
        x=arange(6)
        nfix_obs={'AS':1.95,'ASV':0.53}
        nfix_obs_error={'AS':0.68,'ASV':0.19}
        nfix=columndata['NFIX_TO_SMINN'].sel(time=t_col>(t_col.max()-100)).mean(dim='time')*365*24*3600
        bar(x,nfix,label='N fixation',hatch='//',width=0.4)
        bar(x,columndata['NDEP_TO_SMINN'].sel(time=t_col>(t_col.max()-100)).mean(dim='time')*365*24*3600,bottom=nfix,label='N deposition',hatch='//',width=0.4)
        bar(x+0.4,[nfix_obs.get(e,nan) for e in landscape_ecotypes], yerr=[nfix_obs_error.get(e,nan) for e in landscape_ecotypes]  ,width=0.4,label='Obs N fixation')
        title('N deposition and inputs')
        xticks(x,landscape_ecotypes)
        ylabel('N input rate (gN m$^{-2}$ year$^{-1}$)')
        legend()
        
        subplot(212)
        title('Nutrient and moisture limitation')
        x=arange(len(columndata.lndgrid))
        w=0.8/3
        Nlim=((columndata['GPP']*columndata['FPG']).groupby(year).sum(dim='time')/columndata['GPP'].groupby(year).sum(dim='time'))
        Plim=((columndata['GPP']*columndata['FPG_P']).groupby(year).sum(dim='time')/columndata['GPP'].groupby(year).sum(dim='time'))
        moisturelim=((columndata['GPP']*columndata['BTRAN']).groupby(year).sum(dim='time')/columndata['GPP'].groupby(year).sum(dim='time'))
        bar(x,1-Nlim.mean(dim='year'),yerr=Nlim.std(dim='year'),color='g',label='N limitation (FPG)',width=w,edgecolor='k')
        bar(x+w,1-Plim.mean(dim='year'),yerr=Plim.std(dim='year'),color='orange',label='P limitation (FPG_P)',width=w,edgecolor='k')
        bar(x+w*2,1-moisturelim.mean(dim='year'),yerr=moisturelim.std(dim='year'),color='b',label='Moisture limitation (BTRAN)',width=w,edgecolor='k')
        xticks(x,landscape_ecotypes)
        legend()
        
        tight_layout()
        
        figure('N fixation relations(%s)'%dataname);clf()
        npp=arange(0,200,20)
        
        plot(npp,1.8*(1-exp(-0.003*npp)),'k--',label='Default ELM')
        for pft in range(1,len(pft_names)-1):
            plot(npp,params_new['Nfix_NPP_c1'].values[pft]*(1-exp(-params_new['Nfix_NPP_c2'].values[pft]*npp)),c=pft_colors[pft],label=prettify_pft_name(pft_names[pft]))
        xlabel('NPP (gC m$^{-2}$ year$^{-1}$)')
        ylabel('N fixation (gN m$^{-2}$ year$^{-1}$)')
        legend(fontsize='small')

    years=array([xx.year for xx in vegdata_PFTs.time.values])
    last_100_years=vegdata_PFTs.sel(time=((years.max()-years)<100))

    data_to_plot=last_100_years

    if options.PFTdist or options.all:
        
        figure('PFT distributions',figsize=(9,5));clf()
        plot_PFT_distributions()
        tight_layout()



    if options.biomass or options.all: 
        
        barfig=figure('Biomass comparison by ecotype (%s)'%dataname,figsize=(12,6),clear=True);
        if len(vegdata_PFTs.ecotype)>1:
            axs=barfig.subplots(2,3)
        else:
            axs=barfig.subplots(1,1,squeeze=False)
        
        for econum in range(len(vegdata_PFTs.ecotype)):
            names=[]
            x=0.0
            sca(axs.ravel()[econum])

            h=plot_mod_bar_stack(x+0.4,get_var_PFTs('LEAFC',data_to_plot),econum,width=0.4,hatch='//')
            plot_obs_bar_stack(x,meas_leaf_C.add(meas_nonvasc_C,fill_value=0.0) ,econum,width=0.4)
            names.append('Leaf')
            text(x+0.4,h[-1][0].get_y()+h[-1][0].get_height()+50,'Modeled value',rotation=90,va='bottom',ha='center')
            text(x,h[-1][0].get_y()+h[-1][0].get_height()+50,'Measured value',rotation=90,va='bottom',ha='center')
            x+=1

            plot_mod_bar_stack(x+0.4,get_var_PFTs(['LIVESTEMC','DEADSTEMC'],data_to_plot),econum,width=0.4,hatch='//')
            plot_obs_bar_stack(x,meas_stem_C,econum,width=0.4)
            names.append('Stem')
            x+=1

            plot_mod_bar_stack(x+0.4,get_var_PFTs('FROOTC',data_to_plot),econum,width=0.4,hatch='//')
            bar(x,meas_root_C[landscape_ecotypes[econum]],width=0.4,facecolor=[0.9,0.9,0.9],edgecolor='k')
            names.append('Fine root')
            x+=1

            plot_mod_bar_stack(x+0.4,get_var_PFTs(['LIVECROOTC','DEADCROOTC'],data_to_plot),econum,width=0.4,hatch='//')
            # text(x+0.4,100,'? ? ? ? ? ?',fontsize='large',rotation=90,va='bottom',ha='center')
            plot_obs_bar_stack(x,meas_rhizome_C,econum,width=0.4)
            names.append('Rhizome/\nCoarse root')
            x+=1
            
            plot_mod_bar_stack(x+0.4,get_var_PFTs(['STORVEGC'],data_to_plot),econum,width=0.4,hatch='//')
            names.append('Storage')

            ylabel('Biomass (gC m$^{-2}$)')
            xticks(arange(len(names))+0.3,names)
            title('%s biomass'%ecotype_names_list[econum])
            # ylim(0,3000)
            # ax.set_xlim(right=x+0.6)


        handles=[]
        for pftnum in pfts_inuse[1:]:
            name = prettify_pft_name(data_to_plot.PFTnames.values[pftnum])
            if name == 'Dry Graminoid':
                name = 'Graminoid'
            handles.append(Rectangle([0,0],0,0,facecolor=data_to_plot.PFTcolors.values[pftnum],label=name ))
        l=barfig.axes[0].legend(handles=handles,fontsize='small',ncol=2)
        l.set_draggable(True)

        tight_layout()


        barfig_PFT=figure('Biomass comparison by PFT (%s)'%dataname,figsize=(12,6),clear=True);
        plot_pft_biomass_bars(per_area=False,use_pooled_default_PFTs=True)

        barfig_PFT_perarea=figure('Biomass comparison per area of each PFT (%s)'%dataname,figsize=(12,6),clear=True);
        plot_pft_biomass_bars(per_area=True,use_pooled_default_PFTs=True)
        
        f,a=subplots(num='NPP (%s)'%dataname,figsize=(12,6),clear=True,nrows=2,ncols=3)
        for econum in range(len(data_to_plot.ecotype)):
            ax=a.ravel()[econum]
            x=arange(len(pfts_inuse))
            w=0.8/3
            ax.bar(x,get_var_PFTs('NPP',data_to_plot.sel(PFT=pfts_inuse)).sel(ecotype=econum).mean(dim='time').values*3600*24*365,color=data_to_plot.PFTcolors[pfts_inuse].values,width=w,hatch='//')
            ax.bar(x+w,get_var_PFTs('AGNPP',data_to_plot.sel(PFT=pfts_inuse)).sel(ecotype=econum).mean(dim='time').values*3600*24*365,color=data_to_plot.PFTcolors[pfts_inuse].values,width=w,hatch='/')
            ax.set_title('%s NPP'%ecotype_names_list[econum])
            ax.set_xticks(x)
            ax.set_xticklabels([prettify_pft_name(name) for name in data_to_plot.PFTnames.values[pfts_inuse]],rotation=30,ha='right')
        
        # 
        # w=0.8/6/3
        # for n in range(1,11):
        #     subplot(2,5,n)
        #     title(prettify_pft_name(pft_names[n]))
        #     ylabel('Biomass (gC m$^{-2}$)')
        #     xticks(arange(6)/6,landscape_ecotypes)
        # for pftnum in pfts_inuse[1:]:
        #     if 'c3_arctic_grass' in data_to_plot.PFTnames:
        #         n=pft_names.index(default_new_mappings[data_to_plot.PFTnames.values[pftnum]])
        #     else:
        #         n=pftnum
        #     subplot(2,5,n)
        #     bottom=zeros(6)
        #     bottom_obs=zeros(6)
        #     handles=[]
        # 
        #     vals=get_var_PFTs('FROOTC',data_to_plot,weight_area=False).max(dim='time').sel(PFT=pftnum)
        #     h=bar(arange(6)/6,vals,bottom=bottom,width=w,facecolor='brown',label='Fine root',hatch='//')
        #     bottom=bottom+vals.fillna(0)
        #     handles.append(h)
        # 
        #     vals=get_var_PFTs(['LIVECROOTC','DEADCROOTC'],data_to_plot,weight_area=False).mean(dim='time').sel(PFT=pftnum)
        #     h=bar(arange(6)/6,vals,bottom=bottom,width=w,facecolor='orange',label='Coarse root',hatch='//')
        #     bottom=bottom+vals.fillna(0)
        #     handles.append(h)
        # 
        #     vals=get_var_PFTs('STORVEGC',data_to_plot,weight_area=False).mean(dim='time').sel(PFT=pftnum)
        #     h=bar(arange(6)/6,vals,bottom=bottom,width=w,facecolor='yellow',label='Storage',hatch='//')
        #     bottom=bottom+vals.fillna(0)
        #     handles.append(h)
        #     # bar(arange(6)/6+w,meas_rhizome_C)
        # 
        #     vals=get_var_PFTs(['LIVESTEMC','DEADSTEMC'],data_to_plot,weight_area=False).mean(dim='time').sel(PFT=pftnum)
        #     h=bar(arange(6)/6,vals,bottom=bottom,width=w,facecolor='blue',label='Stem',hatch='//')
        #     bottom=bottom+vals.fillna(0)
        #     handles.append(h)
        # 
        #     vals=get_var_PFTs('LEAFC',data_to_plot,weight_area=False).max(dim='time').sel(PFT=pftnum)
        #     h=bar(arange(6)/6,vals,bottom=bottom,width=w,facecolor='green',label='Leaf',hatch='//')
        #     bottom=bottom+vals.fillna(0)
        #     handles.append(h)  
        # 
        # 
        # 
        # for pftname in unique(meas_leaf_C.index.get_level_values('ELM_PFT')) :
        #     if pftname in ['other','mixed']:
        #          continue
        #     ax=barfig_PFT_perarea.axes[pft_names.index(obsdata_PFT_mappings[pftname])-1]
        #     for econum in range(len(landscape_ecotypes)):
        # 
        #         areafrac=PFT_percents[landscape_ecotypes[econum]][obsdata_PFT_mappings[pftname]]/100
        #         if areafrac==0.0:
        #             invareafrac=0.0
        #         else:
        #             invareafrac=1.0/areafrac
        # 
        #         plotIDs=unique(meas_rhizome_C[landscape_ecotypes[econum]].index.get_level_values('PlotID'))
        #         for plotnum,plotID in enumerate(plotIDs):
        #             bottom=0.0
        #             if pftname not in ['moss','lichen']:
        #                 val=meas_root_C[landscape_ecotypes[econum],plotID]
        #                 h_obsroot=ax.bar(econum/6+w+plotnum*w/len(plotIDs),val,bottom=bottom,facecolor='brown',width=w/len(plotIDs),label='Pooled obs fine roots')
        #                 bottom=bottom+val
        #             val=meas_rhizome_C[landscape_ecotypes[econum],plotID].fillna(0.0).get(pftname,0.0)*invareafrac
        #             h_rhizome=ax.bar(econum/6+w+plotnum*w/len(plotIDs),val,bottom=bottom,facecolor='orange',width=w/len(plotIDs),label='Obs rhizome')
        #             bottom=bottom+val
        #             val=meas_stem_C[landscape_ecotypes[econum],plotID].fillna(0.0).get(pftname,0.0)*invareafrac
        #             ax.bar(econum/6+w+plotnum*w/len(plotIDs),val,bottom=bottom,facecolor='blue',width=w/len(plotIDs))
        #             bottom=bottom+val
        #             val=meas_leaf_C[landscape_ecotypes[econum],plotID].fillna(0.0).get(pftname,0.0)*invareafrac
        #             ax.bar(econum/6+w+plotnum*w/len(plotIDs),val,bottom=bottom,facecolor='green',width=w/len(plotIDs))
        #             bottom=bottom+val
        #             val=meas_nonvasc_C[landscape_ecotypes[econum],plotID].fillna(0.0).get(pftname,0.0)*invareafrac
        #             h_nonvasc=ax.bar(econum/6+w+plotnum*w/len(plotIDs),val,bottom=bottom,facecolor='cyan',width=w/len(plotIDs),label='Obs non-vascular')
        #             bottom=bottom+val
        # 
        # 
        # handles.append(h_nonvasc)
        # handles.append(h_obsroot)
        # handles.append(h_rhizome)
        # l=barfig_PFT_perarea.axes[3].legend(fontsize='small',ncol=1,handles=handles)
        # l.set_draggable(True)
        # 
        # tight_layout()
   
    if options.height or options.all:
        
        figure('Height (%s)'%dataname,figsize=(13,7));clf()


        for econum in range(len(data_to_plot.ecotype)):
            names=[]
            x=0.0
            ax=subplot(2,3,econum+1)
            
            heights=get_var_PFTs('HTOP',data_to_plot,weight_area=False).sel(ecotype=econum).max(dim='time')   
            good=~isnan(heights.values)&(heights.values>0)
            good[data_to_plot.PFTnames=='arctic_wet_graminoid']=False # wet graminoids that we are ignoring
            x=arange(len(heights[good]))
            pftnames=[prettify_pft_name(name) for name in array(data_to_plot.PFTnames)[good]]
            pftcolors=data_to_plot.PFTcolors.values
            if 'arctic_dry_graminoid' in pftnames:
                pftnames[pftnames.index('arctic_dry_graminoid')]='Graminoid'
            for num in range(len(x)):
                bar(x[num],heights[good][num],facecolor=array(pftcolors)[good][num],hatch='//',width=0.4)
                if pftnames[num] in ['Deciduous Shrub Tall','Deciduous Shrub Alder']:
                    bar(x[num]+0.4,obs_heights['tall_shrub_height_mean'],facecolor=array(pftcolors)[good][num],width=0.4)
                elif 'Low' in pftnames[num]:
                    bar(x[num]+0.4,obs_heights['low_shrub_height_mean'],facecolor=array(pftcolors)[good][num],width=0.4)
                elif 'Dwarf' in pftnames[num]:
                    bar(x[num]+0.4,obs_heights['dwarf_shrub_height_mean'],facecolor=array(pftcolors)[good][num],width=0.4)
                elif pftnames[num] is 'Forb':
                    bar(x[num]+0.4,obs_heights['forb_height_mean'],facecolor=array(pftcolors)[good][num],width=0.4)
            xticks(x,pftnames,rotation=45,ha='right')
            ylabel('Height (m)')
            title(ecotype_names_list[econum])
        tight_layout()
   
    if options.AG_BG or options.all:
        
        figure('Aboveground-belowground (%s)'%dataname,figsize=(10.8,4.8));clf()

        aboveground_C = (meas_leaf_C+meas_stem_C).groupby(level=('Ecotype','PlotID')).sum() 
        aboveground_NPP = (meas_leaf_NPP+meas_stem_NPP).groupby(level=('Ecotype','PlotID')).sum()
        belowground_C = meas_rhizome_C.groupby(level=('Ecotype','PlotID')).sum() + meas_root_C
        belowground_NPP = meas_rhizome_NPP.groupby(level=('Ecotype','PlotID')).sum() + meas_froot_NPP

        mod_aboveground_NPP = get_var_PFTs('AGNPP',data_to_plot).mean(dim='time').sum(dim='PFT')*3600*24*365
        mod_belowground_NPP = get_var_PFTs('BGNPP',data_to_plot).mean(dim='time').sum(dim='PFT')*3600*24*365
        mod_aboveground_C  = get_var_PFTs(['LEAFC','LIVESTEMC','DEADSTEMC'],data_to_plot).max(dim='time').sum(dim='PFT')
        mod_belowground_C  = get_var_PFTs(['FROOTC','LIVECROOTC','DEADCROOTC'],data_to_plot).max(dim='time').sum(dim='PFT')

        val=100*(belowground_C/(aboveground_C+belowground_C))[landscape_ecotypes]
        bar(arange(len(landscape_ecotypes)),val.mean(level='Ecotype'),yerr=val.std(level='Ecotype'),width=0.2,label='Biomass')
        val=100*(belowground_NPP/(aboveground_NPP+belowground_NPP))[landscape_ecotypes]
        bar(arange(len(landscape_ecotypes))+0.2,val.mean(level='Ecotype'),yerr=val.std(level='Ecotype'),width=0.2,label='NPP')
        bar(arange(len(landscape_ecotypes))+0.4,100*(mod_belowground_C/(mod_aboveground_C+mod_belowground_C)),width=0.2,label='Mod biomass',facecolor='C0',hatch='//')
        bar(arange(len(landscape_ecotypes))+0.6,100*(mod_belowground_NPP/(mod_aboveground_NPP+mod_belowground_NPP)),width=0.2,label='Mod NPP',facecolor='C1',hatch='//')
        names=[ecotype_names[name] for name in landscape_ecotypes]
        names[0]='Dryas-lichen\ndwarf shrub tundra'
        xticks(arange(len(landscape_ecotypes))+0.3, names,rotation=30,ha='right',fontsize='x-large')
        title('Measured belowground percentage',fontsize='x-large')
        ylabel('Belowground fraction (%)',fontsize='x-large')
        ylim(0,100)
        legend(fontsize='large')
        tight_layout()
   
    if options.nonvasc or options.all:
        
        figure('Non-vascular (%s)'%dataname);clf()
        meas_lichen_BM=Koug_meas_biomass['Nonvascular_gperm2'][:,:,'lichen']
        meas_moss_BM=Koug_meas_biomass['Nonvascular_gperm2'][:,:,'bryophyte']
        # meas_other_nonvasc_BM=Koug_meas_biomass['Nonvascular_gperm2'][:,'other']
        meas_vasc_AG_BM=(Koug_meas_biomass['Leaf_gperm2']+Koug_meas_biomass['Stem_gperm2']).sum(level='PlotID')

        w=0.2
        for econum in range(len(landscape_ecotypes)):
            plotIDs=meas_lichen_BM[landscape_ecotypes[econum]].index.get_level_values('PlotID')
            for plotnum,plotID in enumerate(plotIDs):
                bottom=0.0
                y=meas_lichen_BM[landscape_ecotypes[econum]][plotID]
                lichenbar=bar(econum+plotnum*w/len(plotIDs),y,bottom=bottom,width=w/len(plotIDs),label='Lichen',color=pft_colors[pft_names.index('arctic_lichen')],align='edge')
                bottom=bottom+y
                y=meas_moss_BM[landscape_ecotypes[econum]][plotID]
                mossbar=bar(econum+plotnum*w/len(plotIDs),y,bottom=bottom,width=w/len(plotIDs),label='Moss',color=pft_colors[pft_names.index('arctic_bryophyte')],align='edge')
                y=meas_vasc_AG_BM[plotID]
                vascbar=bar(econum+w+plotnum*w/len(plotIDs),y,facecolor='brown',width=w/len(plotIDs),label='Vascular',align='edge')

        # y=meas_lichen_BM.reindex(landscape_ecotypes,fill_value=0.0)
        # bar(arange(len(landscape_ecotypes)),y,width=0.3,label='Lichen',color=pft_colors[pft_names.index('arctic_lichen')],bottom=bottom)
        # bottom=bottom+y
        # y=meas_moss_BM.reindex(landscape_ecotypes,fill_value=0.0)
        # bar(arange(len(landscape_ecotypes)),y,width=0.3,label='Moss',color=pft_colors[pft_names.index('arctic_bryophyte')],bottom=bottom)
        # bottom=bottom+y
        # # y=meas_other_nonvasc_BM.reindex(landscape_ecotypes,fill_value=0.0)
        # # bar(arange(len(landscape_ecotypes)),y,width=0.3,label='Other nonvasc.',color='gray',bottom=bottom)
        # # bottom=bottom+y
        # y=meas_vasc_AG_BM.reindex(landscape_ecotypes,fill_value=0.0)
        # bar(arange(len(landscape_ecotypes))+0.3,y,width=0.3,label='Vascular',color='brown')
        # # bottom=bottom.add(meas_vasc_AG_BM,fill_value=0.0)

        if 'arctic_bryophyte' in data_to_plot.PFTnames:
            mod_moss=get_var_PFTs('TOTVEGC',data_to_plot,True).max(dim='time').sel(PFT=pft_names.index('arctic_bryophyte'))
            mod_lichen=get_var_PFTs('TOTVEGC',data_to_plot,True).max(dim='time').sel(PFT=pft_names.index('arctic_lichen'))
            mod_vascular = get_var_PFTs(['LEAFC','LIVESTEMC','DEADSTEMC'],data_to_plot).max(dim='time').sel(PFT=arange(3,len(pft_names))).sum(dim='PFT')
        else:
            mod_moss=zeros(len(landscape_ecotypes))
            mod_lichen=zeros(len(landscape_ecotypes))
            mod_vascular = get_var_PFTs(['LEAFC','LIVESTEMC','DEADSTEMC'],data_to_plot).max(dim='time').sum(dim='PFT')
            
        modmossbar=bar(arange(len(landscape_ecotypes))+w*2,mod_lichen,label='Mod moss',color=[.3,1,.3],width=w,hatch='//',align='edge')
        modlichenbar=bar(arange(len(landscape_ecotypes))+w*2,mod_moss,bottom=mod_lichen,label='Mod lichen',color='orange',width=w,hatch='//',align='edge')
        modvascbar=bar(arange(len(landscape_ecotypes))+w*3,mod_vascular,label='Mod vasc',color='brown',width=w,hatch='//',align='edge')
        
        title('Aboveground vascular and non-vascular biomass')
        ylabel('Total biomass (g DW m$^{-2}$)')
        l=legend(handles=(mossbar,lichenbar,vascbar,modmossbar,modlichenbar,modvascbar))
        l.set_draggable(True)
        xticks(arange(len(landscape_ecotypes))+w,[ecotype_names[name] for name in landscape_ecotypes] ,rotation=45,ha='right')
        tight_layout()
       
    if options.vcmax or options.all:
        
        f,axs=subplots(nrows=3,ncols=2,num='VCMAX (%s)'%dataname,clear=True)
        for n in range(len(vegdata_PFTs.ecotype)):
            sca(axs.ravel()[n])
            plot_var_PFTs('VCMAX25TOP',vegdata_PFTs.sel(PFT=pfts_inuse),weight_area=False,ecotype_num=n,maxyear=maxyear,minyear=maxyear-1)
            title(ecotype_names_list[n])
            ylabel('Vcmax(25C)\n(umol m$^{-2}$ s$^{-1}$)')
        legend(fontsize='small',ncol=1)
        # tight_layout()
        
        # x=arange(len(pfts_inuse))
        # bar(x,get_var_PFTs('VCMAX25TOP',data_to_plot.sel(PFT=pfts_inuse),weight_area=False).max(dim='time').mean(dim='ecotype').values,color=data_to_plot.PFTcolors[pfts_inuse].values)
        # xticks(x,[prettify_pft_name(name) for name in data_to_plot.PFTnames.values[pfts_inuse]],rotation=30,ha='right')
        # ylabel('VCmax (umolCO2 m$^{-2}$ s$^{-1}$)')
        # title('Maximum VCmax')
        # subplots_adjust(bottom=0.25)
        
        



    if options.sitemap or options.all: 
        import cartopy.crs as ccrs
        import cartopy
        figure('Site location');clf()
        ax=subplot(111,projection=ccrs.Miller(central_longitude=-155))
        ax.set_extent((-175,-138,50,75)) 
        ax.coastlines(resolution='50m')
        ax.add_feature(cartopy.feature.OCEAN.with_scale('50m'))
        ax.add_feature(cartopy.feature.BORDERS.with_scale('50m'))
        ax.scatter(-164.82,65.16,marker='*',s=50,c='g',transform=ccrs.Geodetic())

    if options.rootprofile or options.all:
        fig=figure('Root profile');clf()
        
        depths=columndata['levdcmp']
        dz=zeros(len(depths))
        dz[0]=0.5*(depths[0]+depths[1])  
        dz[1:-1]=0.5*(depths.shift({'levdcmp':-1})-depths.shift({'levdcmp':1}))[1:-1]
        dz[-1]=depths[-1]-depths[-2]
        
        for eco in vegdata_PFTs['ecotype'].values:
            subplot(2,3,eco+1)
            # for pft in vegdata_PFTs['PFT'][1:-1]:
            for pft in [1,3,10]:
                y=(vegdata_PFTs['FROOT_PROF_unweighted'][:,:,eco,pft].mean(dim='time')*dz).cumsum(dim='levdcmp')
                if y.max()==0:
                    y[:]=nan
                plot(y,columndata['levdcmp'],
                        color=vegdata_PFTs.PFTcolors.values[pft],label=prettify_pft_name(vegdata_PFTs.PFTnames.values[pft]))
            ALTmax=columndata['ALTMAX_LASTYEAR'].isel(lndgrid=eco).mean()
            WTdepth=columndata['ZWT'].isel(lndgrid=eco).mean()
            WTdepth_perched=columndata['ZWT_PERCH'].isel(lndgrid=eco).mean()
            xlims=xlim()
            plot(xlims,[ALTmax,ALTmax],'k--',label='Max active layer depth')
            # plot(xlims,[WTdepth,WTdepth],'b--',label='WT depth')
            plot(xlims,[surfdata['aveDTB'].isel(lsmlon=eco),surfdata['aveDTB'].isel(lsmlon=eco)],'--',label='Bedrock depth')
            plot(xlims,[WTdepth_perched,WTdepth_perched],'b:',label='Perched WT depth')
            ylim(10,-0.5)
            xlabel('Cumulative root fraction')
            ylabel('Depth (m)')
            title(ecotype_names[vegdata_PFTs['community_names'].values[eco]])
        legend(fontsize='small',ncol=1,
            labels=['Nonvascular (beta=0.0)','Evergreen (beta=0.964)','Deciduous (beta=0.914)','Max active layer thickness','Bedrock depth','Perched WT depth'])
        
        tight_layout()
        
        figure('Root depth time series');clf()
        for eco in vegdata_PFTs['ecotype'].values:
            subplot(2,3,eco+1)
            t=array([tt.year+(tt.month-1)/12 for tt in columndata.time.values])
            # for pft in vegdata_PFTs['PFT'][1:-1]:
            for pft in [1,3,10]:
                if (vegdata_PFTs['FROOT_PROF_unweighted'][:,:,eco,pft]*dz).isnull().all():
                    y=zeros(len(t))+nan
                else:
                    y=depths[((vegdata_PFTs['FROOT_PROF_unweighted'][:,:,eco,pft]*dz).cumsum(dim='levdcmp')>0.9).argmax(dim='levdcmp')]

                plot(t,y,
                        color=vegdata_PFTs.PFTcolors.values[pft],label=prettify_pft_name(vegdata_PFTs.PFTnames.values[pft]))
            plot(t,columndata['ALTMAX_LASTYEAR'].isel(lndgrid=eco),'k-')
            # plot(t,columndata['ZWT'].isel(lndgrid=eco),'b-')
            plot([t[0],t[-1]],[surfdata['aveDTB'].isel(lsmlon=eco),surfdata['aveDTB'].isel(lsmlon=eco)],'--',label='Bedrock depth')
            plot(t,columndata['ZWT_PERCH'].isel(lndgrid=eco),'b:')
            title(ecotype_names[vegdata_PFTs['community_names'].values[eco]])
            ylabel('90% cumulative root depth (m)')
            ylims=ylim()
            ylim(ylims[1],ylims[0])
            xlim(1980,2010)
        gcf().axes[0].legend(fontsize='small',ncol=1,
            labels=['Nonvascular (beta=0.0)','Evergreen (beta=0.964)','Deciduous (beta=0.914)','Max active layer thickness','Bedrock depth','Perched WT depth'])
        tight_layout()
   

    if options.save_figs != None:
        save_all_figs(options.save_figs)

    show()

    

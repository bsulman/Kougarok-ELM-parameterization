from kougarok_plotting import *

outputdata_dir='../../output_data'


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


vegdata_spinup = read_pftfile(outputdata_dir+'/ELMuserpft_Kougarok_rhizomes-as-storage_h1_20190418.nc',maxyear=400)
columndata_spinup = xarray.open_dataset(outputdata_dir+'/ELMuserpft_Kougarok_rhizomes-as-storage_h0_20190418.nc',autoclose=True)

vegdata_oldparams_hist = read_pftfile(outputdata_dir+'/ELMuserpft_Kougarok_ICB20TRCNPRDCTCBC.clm2.h_PFTs.nc')
spinup_oldparams=read_pftfile('../../output_data/ELMuserpft_adspinuptest_Kougarok_ICB1850CNPRDCTCBC.h1.nc')

Koug_meas_biomass.rename(index={'NAMC':'DL','TTWBT':'ASV'},inplace=True)
Koug_meas_chem.rename(index={'NAMC':'DL','TTWBT':'ASV'},inplace=True)

meas_leaf_C=(Koug_meas_biomass['LeafBiomass_gperm2']*Koug_meas_chem['LeafC_percent']/100)
meas_nonvasc_C=(Koug_meas_biomass['NonvascularBiomass_gperm2']*0.5)
meas_rhizome_C=(Koug_meas_biomass['RhizomeBiomass_gperm2']*Koug_meas_chem['RhizomeC_percent']/100).groupby('Ecotype').sum()
meas_root_C=Koug_meas_biomass['FineRootBiomass_gperm2'][:,'mixed']*(Koug_meas_chem['FineRootC_percent']/100)
meas_stem_C=(Koug_meas_biomass['StemBiomass_gperm2']*Koug_meas_chem['StemC_percent']/100)
meas_rhizome_C=(Koug_meas_biomass['RhizomeBiomass_gperm2']*Koug_meas_chem['RhizomeC_percent']/100)
meas_NPP_C=(Koug_meas_biomass['LeafNPP_gperm2peryr']*Koug_meas_chem['LeafC_percent']/100)+\
    (Koug_meas_biomass['StemNPP_gperm2peryr']*Koug_meas_chem['StemC_percent']/100)#+\

meas_froot_NPP =(Koug_meas_biomass['FineRootNPP_gperm2peryr']*Koug_meas_chem['FineRootC_percent']/100).groupby('Ecotype').sum()


figure('PFT distributions');clf()
subplot(121)
names=[]
bottom=zeros(len(landscape_ecotypes))
for pftnum in range(len(pft_names_default[:17])):
    pft_pcts=PFT_percents_default.loc[pft_names_default[pftnum]]
    if (pft_pcts==0).all():
        continue
    bar(arange(len(landscape_ecotypes)),pft_pcts,bottom=bottom)
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



def plot_mod_bar_stack(x,dat,econum,do_legend=False,**kwargs):
    handles=[]
    bottom=0.0
    for pftnum in range(len(pft_names)):
        if pftnum==10:
            val=dat.sel(PFT=10,ecotype=econum).max(dim='time')+dat.sel(PFT=11,ecotype=econum).max(dim='time')
        elif pftnum==11:
            continue
        else:
            val=dat.sel(PFT=pftnum,ecotype=econum).max(dim='time')
        if ~isnan(val):
            if pft_names[pftnum].startswith('arctic'):
                name=pft_names[pftnum][len('arctic_'):]
            else:
                name=pft_names[pftnum]
            if name=='dry_graminoid':
                name='graminoid'
            handles.append(bar(x,val,bottom=bottom,facecolor=pft_colors[pftnum],label=name,**kwargs))
            bottom=bottom+val
    return handles
    


def plot_obs_bar_stack(x,obsdata,ecotype_num,**kwargs):
    names=[]
    bottom=0.0
    for pft in obsdata[landscape_ecotypes[ecotype_num]].index:
        val=obsdata[(landscape_ecotypes[ecotype_num],pft)]

        if ~isnan(val):
            bar(x,val,bottom=bottom,facecolor=pft_colors[(pft_names.index(obsdata_PFT_mappings[pft]))],**kwargs)
            bottom=bottom+val

def prettify_pft_name(name):
    if name.startswith('arctic'):
        pretty_name=name[len('arctic_'):]
    else:
        pretty_name=name
    pretty_name = ' '.join(pretty_name.split('_')).title()
    return pretty_name

barfig=figure('Biomass comparison',figsize=(15,8));clf()
for econum in range(6):
    names=[]
    x=0.0
    ax=subplot(2,3,econum+1)

    h=plot_mod_bar_stack(x+0.4,get_var_PFTs('LEAFC',vegdata_oldparams_hist),econum,width=0.4,hatch='//')
    plot_obs_bar_stack(x,meas_leaf_C.add(meas_nonvasc_C,fill_value=0.0) ,econum,width=0.4)
    names.append('Leaf')
    text(x+0.4,h[-1][0].get_y()+h[-1][0].get_height()+50,'Modeled value',rotation=90,va='bottom',ha='center')
    text(x,h[-1][0].get_y()+h[-1][0].get_height()+50,'Measured value',rotation=90,va='bottom',ha='center')
    x+=1

    plot_mod_bar_stack(x+0.4,get_var_PFTs(['LIVESTEMC','DEADSTEMC'],vegdata_oldparams_hist),econum,width=0.4,hatch='//')
    plot_obs_bar_stack(x,meas_stem_C,econum,width=0.4)
    names.append('Stem')
    x+=1

    plot_mod_bar_stack(x+0.4,get_var_PFTs('FROOTC',vegdata_oldparams_hist),econum,width=0.4,hatch='//')
    bar(x,meas_root_C[landscape_ecotypes[econum]],width=0.4,facecolor=[0.9,0.9,0.9],edgecolor='k')
    names.append('Fine root')
    x+=1

    plot_mod_bar_stack(x+0.4,get_var_PFTs(['LIVECROOTC','DEADCROOTC'],vegdata_oldparams_hist),econum,width=0.4,hatch='//')
    # text(x+0.4,100,'? ? ? ? ? ?',fontsize='large',rotation=90,va='bottom',ha='center')
    plot_obs_bar_stack(x,meas_rhizome_C,econum,width=0.4)
    names.append('Rhizome/\nCoarse root')
    x+=1

    ylabel('Biomass (gC m$^{-2}$)')
    xticks(arange(len(names))+0.3,names)
    title('%s biomass'%ecotype_names_list[econum])
    ylim(0,3000)
    # ax.set_xlim(right=x+0.6)

handles=[]
for pftnum in range(1,len(pft_names)-1):
    name = prettify_pft_name(pft_names[pftnum])
    if name == 'Dry Graminoid':
        name = 'Graminoid'
    handles.append(Rectangle([0,0],0,0,facecolor=pft_colors[pftnum],label=name ))
l=barfig.axes[1].legend(handles=handles,fontsize='small',ncol=2)
l.set_draggable(True)

tight_layout()





figure('Height');clf()
heights=get_var_PFTs('HTOP',vegdata_oldparams_hist,weight_area=False).sel(ecotype=2).max(dim='time')   
good=~isnan(heights.values)&(heights.values>0)
good[11]=False # wet graminoids that we are ignoring
x=arange(len(heights[good]))
pftnames=[prettify_pft_name(name) for name in array(pft_names)[good]]
pftnames[-1]='Graminoid'
for num in range(len(x)):
    bar(x[num],heights[good][num],facecolor=array(pft_colors)[good][num])
xticks(x,pftnames,rotation=45,ha='right')
ylabel('Height (m)')
title('Modeled vegetation height')
tight_layout()

figure('Aboveground-belowground',figsize=(10.8,4.8));clf()
meas_leaf_C=(Koug_meas_biomass['LeafBiomass_gperm2']*Koug_meas_chem['LeafC_percent']/100)
meas_root_C=(Koug_meas_biomass['FineRootBiomass_gperm2'][:,'mixed']*(Koug_meas_chem['FineRootC_percent']/100))[:,'mixed']
meas_rhizome_C=(Koug_meas_biomass['RhizomeBiomass_gperm2']*Koug_meas_chem['RhizomeC_percent']/100)
meas_rhizome_NPP=(Koug_meas_biomass['RhizomeNPP_gperm2peryr']*Koug_meas_chem['RhizomeC_percent']/100)
meas_leaf_NPP   =(Koug_meas_biomass['LeafNPP_gperm2peryr']*Koug_meas_chem['LeafC_percent']/100)
meas_stem_NPP   =(Koug_meas_biomass['StemNPP_gperm2peryr']*Koug_meas_chem['StemC_percent']/100)
meas_root_NPP   =(Koug_meas_biomass['FineRootNPP_gperm2peryr']*Koug_meas_chem['FineRootC_percent']/100)[:,'mixed']

aboveground_C = (meas_leaf_C+meas_stem_C).sum(level='Ecotype')
aboveground_NPP = (meas_leaf_NPP+meas_stem_NPP).sum(level='Ecotype')
belowground_C = meas_rhizome_C.sum(level='Ecotype') + meas_root_C
belowground_NPP = meas_rhizome_NPP.sum(level='Ecotype') + meas_root_NPP

mod_aboveground_NPP = get_var_PFTs('AGNPP',vegdata_oldparams_hist).mean(dim='time').sum(dim='PFT')*3600*24*365
mod_belowground_NPP = get_var_PFTs('BGNPP',vegdata_oldparams_hist).mean(dim='time').sum(dim='PFT')*3600*24*365
mod_aboveground_C  = get_var_PFTs(['LEAFC','LIVESTEMC','DEADSTEMC'],vegdata_oldparams_hist).max(dim='time').sum(dim='PFT')
mod_belowground_C  = get_var_PFTs(['FROOTC','LIVECROOTC','DEADCROOTC'],vegdata_oldparams_hist).max(dim='time').sum(dim='PFT')

bar(arange(len(landscape_ecotypes)),100*(belowground_C/(aboveground_C+belowground_C))[landscape_ecotypes],width=0.2,label='Biomass')
bar(arange(len(landscape_ecotypes))+0.2,100*(belowground_NPP/(aboveground_NPP+belowground_NPP))[landscape_ecotypes],width=0.2,label='NPP')
bar(arange(len(landscape_ecotypes))+0.4,100*(mod_belowground_C/(mod_aboveground_C+mod_belowground_C)),width=0.2,label='Mod biomass',facecolor='C0',edgecolor='k')
bar(arange(len(landscape_ecotypes))+0.6,100*(mod_belowground_NPP/(mod_aboveground_NPP+mod_belowground_NPP)),width=0.2,label='Mod NPP',facecolor='C1',edgecolor='k')
names=[ecotype_names[name] for name in landscape_ecotypes]
names[0]='Dryas-lichen\ndwarf shrub tundra'
xticks(arange(len(landscape_ecotypes))+0.3, names,rotation=30,ha='right',fontsize='x-large')
title('Measured belowground percentage',fontsize='x-large')
ylabel('Belowground fraction (%)',fontsize='x-large')
ylim(0,100)
legend(fontsize='large')
tight_layout()


figure('Non-vascular');clf()
meas_lichen_BM=Koug_meas_biomass['NonvascularBiomass_gperm2'][:,'lichen']
meas_moss_BM=Koug_meas_biomass['NonvascularBiomass_gperm2'][:,'moss']
meas_other_nonvasc_BM=Koug_meas_biomass['NonvascularBiomass_gperm2'][:,'other']
meas_vasc_AG_BM=(Koug_meas_biomass['LeafBiomass_gperm2']+Koug_meas_biomass['StemBiomass_gperm2']).sum(level='Ecotype')

bottom=meas_lichen_BM.reindex(landscape_ecotypes,fill_value=0.0)*0

y=meas_lichen_BM.reindex(landscape_ecotypes,fill_value=0.0)
bar(arange(len(landscape_ecotypes)),y,width=0.3,label='Lichen',color=pft_colors[pft_names.index('arctic_lichen')],bottom=bottom)
bottom=bottom+y
y=meas_moss_BM.reindex(landscape_ecotypes,fill_value=0.0)
bar(arange(len(landscape_ecotypes)),y,width=0.3,label='Moss',color=pft_colors[pft_names.index('arctic_bryophyte')],bottom=bottom)
bottom=bottom+y
y=meas_other_nonvasc_BM.reindex(landscape_ecotypes,fill_value=0.0)
bar(arange(len(landscape_ecotypes)),y,width=0.3,label='Other nonvasc.',color='gray',bottom=bottom)
bottom=bottom+y
y=meas_vasc_AG_BM.reindex(landscape_ecotypes,fill_value=0.0)
bar(arange(len(landscape_ecotypes))+0.3,y,width=0.3,label='Vascular',color='brown')
# bottom=bottom.add(meas_vasc_AG_BM,fill_value=0.0)

mod_moss=get_var_PFTs('TOTVEGC',vegdata_oldparams_hist,0).max(dim='time').sel(PFT=pft_names.index('arctic_bryophyte'))
mod_lichen=get_var_PFTs('TOTVEGC',vegdata_oldparams_hist,0).max(dim='time').sel(PFT=pft_names.index('arctic_lichen'))
bar(arange(len(landscape_ecotypes))+0.3+0.3,mod_lichen,label='Mod moss',color=[.3,1,.3],width=0.3)
bar(arange(len(landscape_ecotypes))+0.3+0.3,mod_moss,bottom=mod_lichen,label='Mod lichen',color='orange',width=0.3)
title('Aboveground vascular and non-vascular biomass')
ylabel('Total biomass (g DW m$^{-2}$)')
l=legend()
l.set_draggable(True)
xticks(arange(len(landscape_ecotypes))+0.3,[ecotype_names[name] for name in landscape_ecotypes] ,rotation=45,ha='right')
tight_layout()

import cartopy.crs as ccrs

figure('Site location');clf()
ax=subplot(111,projection=ccrs.Miller(central_longitude=-155))
ax.set_extent((-175,-138,50,75)) 
ax.coastlines(resolution='50m')
ax.add_feature(cartopy.feature.OCEAN.with_scale('50m'))
ax.add_feature(cartopy.feature.BORDERS.with_scale('50m'))
ax.scatter(-164.82,65.16,marker='*',s=50,c='g',transform=ccrs.Geodetic())


figure('Evergreen vs decid');clf()
econum=1
leaf_decid=get_var_PFTs('LEAFC',data).isel(PFT=10,ecotype=econum)+get_var_PFTs('LEAFC',data).isel(PFT=11,ecotype=econum)
froot_decid=get_var_PFTs('FROOTC',data).isel(PFT=10,ecotype=econum)+get_var_PFTs('FROOTC',data).isel(PFT=11,ecotype=econum)
storage_decid=get_var_PFTs('STORVEGC',data).isel(PFT=10,ecotype=econum)+get_var_PFTs('STORVEGC',data).isel(PFT=11,ecotype=econum)

t_decid=array([tt.year + (tt.dayofyr-1)/365 for tt in leaf_decid['time'].data])
leaf_eg=get_var_PFTs('LEAFC',vegdata_spinup).isel(PFT=10,ecotype=econum)
froot_eg=get_var_PFTs('FROOTC',vegdata_spinup).isel(PFT=10,ecotype=econum)
storage_eg=get_var_PFTs('STORVEGC',vegdata_spinup).isel(PFT=10,ecotype=econum)
t_eg=array([tt.year + (tt.dayofyr-1)/365 for tt in leaf_eg['time'].data])
subplot(121)
plot(t_decid,leaf_decid,label='Leaf',c='g')
plot(t_decid,froot_decid,label='Fine-root',c='C0')
plot(t_decid,storage_decid,label='Storage',c='C1')
xlim(80,100)
title('Deciduous graminoid',fontsize='x-large')
xlabel('Time (years)',fontsize='x-large')
ylabel('Biomass (gC m$^{-2}$)',fontsize='x-large')
# ylim(-0.1,12)
legend(fontsize='x-large')
gca().set_ylim(bottom=0)
xticks(fontsize='large')
yticks(fontsize='large')

subplot(122)
plot(t_eg,leaf_eg,label='Leaf',c='g')
plot(t_eg,froot_eg,label='Fine-root',c='C0')
plot(t_eg,storage_eg,label='Storage',c='C1')
xlim(80,100)
title('Evergreen graminoid',fontsize='x-large')
xlabel('Time (years)',fontsize='x-large')
ylabel('Biomass (gC m$^{-2}$)',fontsize='x-large')
gca().set_ylim(bottom=0)
xticks(fontsize='large')
yticks(fontsize='large')

tight_layout()

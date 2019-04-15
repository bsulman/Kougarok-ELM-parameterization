from kougarok_plotting import *

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

outputdata_dir='../../output_data'
vegdata_PFTs_oldparams=read_pftfile(outputdata_dir+'/ELMuserpft_adspinuptest_Kougarok_ICB1850CNPRDCTCBC.h1.nc')
vegdata_PFTs_newparams=read_pftfile(outputdata_dir+'/ELMuserpft_adspinuptest_noPlim_newparams.h1_20190308.nc',maxyear=150)

maxyear=50
meas_leaf_C=(Koug_meas_biomass['LeafBiomass_gperm2']*Koug_meas_chem['LeafC_percent']/100).copy()
meas_nonvasc_C=(Koug_meas_biomass['NonvascularBiomass_gperm2']*0.5)
# Should this include nonvascular biomass?
for ecotype in landscape_ecotypes:
    meas_leaf_C[ecotype,'moss']=meas_nonvasc_C.loc[ecotype,'moss']
    meas_leaf_C[ecotype,'lichen']=meas_nonvasc_C.loc[ecotype,'lichen']
heightfig=plot_pair('HTOP',vegdata_PFTs_oldparams,vegdata_PFTs_newparams,weight_area=False,maxyear=maxyear)
leaffig=plot_pair('LEAFC',vegdata_PFTs_oldparams,vegdata_PFTs_newparams,obsdata=meas_leaf_C,maxyear=maxyear)
meas_root_C=Koug_meas_biomass['FineRootBiomass_gperm2'][:,'mixed']*(Koug_meas_chem['FineRootC_percent']/100)
frootfig=plot_pair('FROOTC',vegdata_PFTs_oldparams,vegdata_PFTs_newparams,obsdata=meas_root_C,maxyear=maxyear,plotsum=True)



figure('Temperature and root respiration cumulative');clf()
Tsoil10cm=xarray.open_dataset(outputdata_dir+'/ELMuserpft_Kougarok_ICB1850CNPRDCTCBC_clm2_h_20190129.nc',autoclose=True)['TSOI_10CM']
t2=array([tt.year + (tt.month-.5)/12 for tt in Tsoil10cm.time.data])
plot(t2,Tsoil10cm.isel(lndgrid=1)-273.15,'b-')
plot([0,maxyear],[0.0,0.0],'k--')
ylabel('Soil temperature (C)')

ax2=twinx()
t=array([tt.year + (tt.dayofyr-1)/365 for tt in vegdata_PFTs_newparams['time'].data])
ecotype_num=0
startyear=15
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

# vegdata_PFTs_dormancy=read_pftfile(outputdata_dir+'/ELMuserpft_adspinuptest_newparams_Kougarok_ICB1850CNPRDCTCBC.h1_dormantMR_2019-02-26.nc')
vegdata_PFTs_dormancy=read_pftfile(outputdata_dir+'/ELMuserpft_adspinuptest_noPlim_newparams.h1_20190308.nc')

figure('Temperature and root respiration cumulative dormancy');clf()
# Tsoil10cm=xarray.open_dataset(outputdata_dir+'/ELMuserpft_Kougarok_ICB1850CNPRDCTCBC_clm2_h_20190129.nc',autoclose=True)['TSOI_10CM']
t2=array([tt.year + (tt.month-.5)/12 for tt in Tsoil10cm.time.data])
plot(t2,Tsoil10cm.isel(lndgrid=1)-273.15,'b-')
plot([0,maxyear],[0.0,0.0],'k--')
ylabel('Soil temperature (C)')

ax2=twinx()
t=array([tt.year + (tt.dayofyr-1)/365 for tt in vegdata_PFTs_dormancy['time'].data])
ecotype_num=0
startyear=15
endyear=25
for yr in range(startyear,endyear):
    growingseason_start=nonzero((t>=yr)&(vegdata_PFTs_dormancy.sel(PFT=3,ecotype=ecotype_num)['GPP_unweighted'].values>0))[0][0]
    growingseason_start_nextyear=nonzero((t>=yr+1)&(vegdata_PFTs_dormancy.sel(PFT=3,ecotype=ecotype_num)['GPP_unweighted'].values>0))[0][0]
    xx=arange(growingseason_start,growingseason_start_nextyear)
    plot_var_PFTs('MR',vegdata_PFTs_dormancy.isel(time=xx),ecotype_num,cumulative=True,ax=ax2,modfactor=3600*24,units='gC/m2',minyear=yr,ls='--')
    plot_var_PFTs('FROOT_MR',vegdata_PFTs_dormancy.isel(time=xx),ecotype_num,cumulative=True,ax=ax2,modfactor=3600*24,units='gC/m2',minyear=yr,ls=':')
    plot_var_PFTs('GPP',vegdata_PFTs_dormancy.isel(time=xx),ecotype_num,cumulative=True,ax=ax2,modfactor=3600*24,units='gC/m2')
xlim(startyear,endyear)

tight_layout()


figure('XSMR pool');clf()
plot_var_PFTs('XSMRPOOL',vegdata_PFTs_dormancy,ecotype_num,maxyear=30,longname='Maint. Resp imbalance pool')
ylim(-320,30)
# title('Maint. Resp imbalance pool')

figdir='../Figures/2019-02-27'

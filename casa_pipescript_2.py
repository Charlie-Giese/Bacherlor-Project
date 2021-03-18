__rethrow_casa_exceptions = True
context = h_init()
context.set_state('ProjectSummary', 'observatory', 'Karl G. Jansky Very Large Array')
context.set_state('ProjectSummary', 'telescope', 'EVLA')

data_file='19B-105.sb37462090.eb37534141.58806.000772141204'

try:


    hifv_importdata(vis=[data_file], createmms='automatic', 
    				asis='Receiver CalAtmosphere', ocorr_mode='co', 
    				nocopy=False, overwrite=False)


    hifv_hanning(pipelinemode="automatic")


    hifv_flagdata(tbuff=0.0, flagbackup=False, scan=True, fracspw=0.05, intents='*POINTING*,*FOCUS*,*ATMOSPHERE*,*SIDEBAND_RATIO*, *UNKNOWN*, *SYSTEM_CONFIGURATION*, *UNSPECIFIED#UNSPECIFIED*', clip=True, baseband=True, shadow=True, quack=True, edgespw=True, autocorr=True, hm_tbuff='1.5int', template=True, online=True)
    flagdata(vis=data_file+'.ms', antenna='ea16', flagbackup=False)   #add command for each antenna I want to flag
    flagdata(vis=data_file+'.ms', antenna='ea19', flagbackup=False)
    flagdata(vis=data_file+'.ms', antenna='ea8', flagbackup=False)
    flagdata(vis=data_file+'.ms', spw='65', flagbackup=False)
    flagdata(vis=data_file+'.ms', spw='67', flagbackup=False)
    flagdata(vis=data_file+'.ms', spw='79', flagbackup=False)
    flagdata(vis=data_file+'.ms', antenna='ea04', spw='32,47', flagbackup=False)
    flagdata(vis=data_file+'.ms', antenna='ea05', spw='17,20', flagbackup=False)
    flagdata(vis=data_file+'.ms', antenna='ea06', spw='47', flagbackup=False)
    flagdata(vis=data_file+'.ms', antenna='ea07', spw='17,20', flagbackup=False)
    flagdata(vis=data_file+'.ms', antenna='ea13', spw='17,32,41,47,61', flagbackup=False)
    flagdata(vis=data_file+'.ms', antenna='ea14', spw='47,77', flagbackup=False)
    flagdata(vis=data_file+'.ms', antenna='ea22', spw='32,48,61', flagbackup=False)
    flagdata(vis=data_file+'.ms', antenna='ea23', spw='61', flagbackup=False)
    flagdata(vis=data_file+'.ms', antenna='ea27', spw='17,32,41', flagbackup=False)

    hifv_vlasetjy(fluxdensity=-1, scalebychan=True, spix=0, reffreq='1GHz')
    hifv_priorcals(tecmaps=False)
    hifv_testBPdcals(weakbp=False, refantignore='ea12')
    hifv_flagbaddef(doflagundernspwlimit=True)
    hifv_checkflag(pipelinemode="automatic")
    hifv_semiFinalBPdcals(weakbp=False, refantignore='ea12')
    hifv_checkflag(checkflagmode='semi')
    hifv_semiFinalBPdcals(weakbp=False, refantignore='ea12')
    hifv_solint(pipelinemode="automatic", refantignore='ea12')
    hifv_fluxboot2(fitorder=-1, refantignore='ea12')
    hifv_finalcals(weakbp=False, refantignore='ea12')
    hifv_applycals(flagdetailedsum=True, gainmap=False, flagbackup=True, flagsum=True)
    hifv_targetflag(intents='*CALIBRATE*,*TARGET*')
    hifv_statwt(datacolumn='corrected')
    hifv_plotsummary(pipelinemode="automatic")
    #hif_makeimlist(nchan=-1, calcsb=False, intent='PHASE,BANDPASS', robust=-999.0, parallel='automatic', per_eb=False, calmaxpix=300, specmode='cont', clearlist=True)
    #hif_makeimages(tlimit=2.0, hm_perchanweightdensity=False, hm_npixels=0, hm_dogrowprune=True, hm_negativethreshold=-999.0, calcsb=False, hm_noisethreshold=-999.0, hm_fastnoise=True, hm_masking='none', hm_minpercentchange=-999.0, parallel='automatic', masklimit=4, hm_nsigma=0.0, target_list={}, hm_minbeamfrac=-999.0, hm_lownoisethreshold=-999.0, hm_growiterations=-999, overwrite_on_export=True, cleancontranges=False, hm_sidelobethreshold=-999.0)
    hifv_exportdata(gainmap=False, exportmses=False, exportcalprods=False)
finally:
    h_save()

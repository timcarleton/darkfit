import astroquery
from astroquery.mast import Observations
from astropy import units, coordinates
import numpy as np

#narrowfilters=['F126N','F128N','F130N','F132N']
#widefilter='F125W'
#narrowfilters=['F164N','F167N']
#widefilter='F160W'
narrowfilters=['F160W']
widefilter='F125W'

#rather than using query, use existing image lists

for narrowfilter in narrowfilters:
    obsnarrow=Observations.query_criteria(obs_collection='HST',filters=[narrowfilter],t_exptime=[200,20000],intentType='science',calib_level=3)
    obsnarrowf=Observations.filter_products(Observations.get_product_list(obsnarrow),productSubGroupDescription='FLT',type='S')

    narrowcoords=coordinates.SkyCoord(obsnarrow['s_ra']*units.deg,obsnarrow['s_dec']*units.deg).galactic
    wextragalactic=np.where(abs(narrowcoords.b)>20*units.deg)[0]
    wextragalactic=np.array([i for i in wextragalactic if (obsnarrow['obsid'][i] in obsnarrowf['parent_obsid'])])
    #wextragalactic=np.array([i for i in wextragalactic if ('flt.fits' in obsnarrow[i]['dataURL'])])

    matchednarrow=[]
    matchedwide=[]
    
    print(wextragalactic)
    for i in wextragalactic:
        obswide=Observations.query_criteria(obs_collection='HST',filters=[widefilter],proposal_id=[obsnarrow[i]['proposal_id']],s_ra=[obsnarrow[i]['s_ra']-.016*2,obsnarrow[i]['s_ra']+.016*2],s_dec=[obsnarrow[i]['s_dec']-.016*2,obsnarrow[i]['s_dec']+.016*2])
        print(i,len(obswide))
        if len(obswide)>0:
            for j in obswide['obs_id']:
                if len(j)==9:
                    print(obsnarrow[i]['obs_id'],obswide['obs_id'])
                    matchedwide.append(j)
            matchednarrow.append(obsnarrow[i]['obs_id'])

    np.savetxt('narrowids_wide'+narrowfilter+'.txt',np.array([matchednarrow]).T,'%s')
    np.savetxt('wideids_wide'+narrowfilter+'_'+widefilter+'.txt',np.array([matchedwide]).T,fmt='%s')

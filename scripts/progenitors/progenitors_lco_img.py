import util
import sheetproc
import sys
import datetime
import numpy as np
from lcogt import lcogt
from astropy.time import Time, TimeDelta
from astropy.table import unique, vstack, Table
from astropy.coordinates import SkyCoord
from astropy import units as u

band_priority=['r','rp','r-ZTF','R','orange','V','i','ip','I','g','gp','g-ZTF',
    'cyan','B']
lco = lcogt('/home/ckilpatrick/scripts/shibboleth')
sdate = Time(datetime.datetime.now()) - TimeDelta(28, format='jd')
rg_img = lco.get_requestgroups(propid=lco.proposals,
    itype='2M0-SCICAM-SPECTRAL,2M0-SCICAM-MUSCAT', sdate=sdate)

def main(redo=False, download_yse=True):

    # Build out progenitor data from google sheet and add metadata from files
    sndata = sheetproc.download_progenitor_data(util.params['SHEET'])
    snphot = util.get_yse_target_photometry()

    use_targets = None
    now = Time(datetime.datetime.now()).mjd

    # Iterate through data by spectroscopic type - Ia, IIb, IIn, II, Ib/c, Other
    all_keys = list(sndata.keys())
    for type_key in all_keys:
        # We're not doing SN Ia follow up here unless in specific cases
        if 'Other' in type_key: continue
        type_data = sndata[type_key]

        # Define initial target sample
        mask=now-Time(type_data['Discovery Date'].data).mjd<800
        if not use_targets:
            use_targets = type_data[mask]
        else:
            use_targets = vstack([use_targets, type_data[mask]])

    for target in use_targets:
        mask = snphot['name']==target['Name']
        photdata = snphot[mask]

        # Want to assess whether target is detectable based on optical data
        # BVgri
        mask = np.array([r['filter'] in band_priority for r in photdata])
        photdata = photdata[mask]

        mask = now-Time(photdata['obsdate']).mjd < 14
        photdata = photdata[mask]

        if not photdata or len(photdata)==0: continue

        # Determine most recent brightness in following band priority
        mag=0.0
        for bp in band_priority:
            if len(photdata[photdata['filter']==bp])>0:
                subphot = photdata[photdata['filter']==bp]
                subphot.sort('obsdate')
                mag=subphot['mag'][-1]
                magerr=subphot['magerr'][-1]
                filt=bp
                break

        if mag>23: continue

        needs_obs = True
        cadence = 7
        targcoord = SkyCoord(target['RA'], target['Dec'], unit=(u.hour,u.deg))

        # Get LCO data to determine if we need a new observation
        for request in rg_img:
            # If PENDING observation, don't need a new one
            if request['state']=='PENDING':
                for config in request['requests'][0]['configurations']:
                    print(config['target'])
                    if ('ra' not in config['target'].keys() or
                        'dec' not in config['target'].keys()):
                        continue
                    print(config['target']['ra'],config['target']['dec'])
                    coord = SkyCoord(config['target']['ra'],
                        config['target']['dec'], unit='deg')
                    if coord.separation(targcoord).degree < 0.3:
                        needs_obs = False

            # Use start time as a proxy for observation time and cut on observations
            # that are older than the now - cadence
            time = Time(request['requests'][0]['modified']).mjd
            now = Time(datetime.datetime.now()).mjd
            if now - time > cadence:
                continue

            # For records more recent than cadence, check COMPLETED and PENDING
            # observations.
            if request['state']=='COMPLETED':
                for config in request['requests'][0]['configurations']:
                    coord = SkyCoord(config['target']['ra'],
                        config['target']['dec'], unit='deg')
                    # If the separation between the recent or pending observation
                    # and the target coordinates is smaller than threshold, then
                    # record that this target does not need to be re-observed
                    if coord.separation(targcoord).degree < 0.3:
                        needs_obs = False

        if needs_obs:
            ra = targcoord.ra.degree
            dec = targcoord.dec.degree

            propidx=0
            try_next_proposal = True

            while try_next_proposal:

                if dec > -30:
                    print('Trying MuSCAT3')
                    response = lco.make_obs_request(target['Name'], ra, dec,
                        mag, propidx=propidx, strategy = 'photometry-muscat')
                else:
                    print('Trying Spectral')
                    response = lco.make_obs_request(target['Name'], ra, dec,
                        mag, propidx=propidx, strategy = 'photometry-spectral')

                if response:
                    if 'non_field_errors' in response.keys():
                        resp = response['non_field_errors'][0]
                        if ('time' in resp and 'allocated' in resp):
                            print(response['non_field_errors'])
                            propidx += 1
                            continue

                try_next_proposal=False

            if response and 'requests' in response.keys():
                if 'non_field_errors' in response['requests'][0].keys():
                    message = '{target} cannot be scheduled due to availability.'
                    print(response['requests'][0])
                    print(message.format(target=target['Name']))
                else:
                    message = 'Successfully scheduled {target}'
                    print(message.format(target=target['Name']))
            else:
                message = '{target} exposures not possible with current settings'
                print(message.format(target=target['Name']))
                print(response)
        else:
            message = '{target} does not need a new observation.'
            print(message.format(target=target['Name']))

main(redo=False)


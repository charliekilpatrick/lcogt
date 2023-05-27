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
rg_spec = lco.get_requestgroups(propid=lco.proposals,
    itype='2M0-FLOYDS-SCICAM', sdate=sdate)
rg_img = lco.get_requestgroups(propid=lco.proposals,
    itype='2M0-SCICAM-SPECTRAL,2M0-SCICAM-MUSCAT', sdate=sdate)

def schedule_target(target, snphot, now):
        mask = snphot['name']==target['Name']
        photdata = snphot[mask]

        # Want to assess whether target is detectable based on optical data
        # BVgri
        mask = np.array([r['filter'] in band_priority for r in photdata])
        photdata = photdata[mask]

        mask = now-Time(photdata['obsdate']).mjd < 14
        photdata = photdata[mask]

        if not photdata or len(photdata)==0:
            print('No phot data for {0}'.format(target['Name']))
            return(0)

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

        if mag>18.5: return(0)

        needs_obs = True
        cadence = 1
        targcoord = SkyCoord(target['RA'], target['Dec'], unit=(u.hour,u.deg))
        lco.params['strategy']['default']['ipp']=1.2
        lco.params['strategy']['default']['window']=1.0

        # Get LCO data to determine if we need a new observation
        for request in rg_spec:
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

                response = lco.make_obs_request(target['Name'], ra, dec,
                    mag, propidx=propidx, strategy = 'spectroscopy')

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
        schedule_target(target, snphot, now)

    # Adding other high-priority targets
    targets = [{'Name': '2021oat', 'RA': 195.034453529,'Dec':28.1703005032}]
    for target in targets:
        schedule_target(target, snphot, now)

main(redo=False)


#!/usr/bin/env python3
from lcogt import lcogt
from datetime import tzinfo, timedelta, datetime
from astropy import units as u
from astropy.table import Table,unique
from astropy.coordinates import SkyCoord
from astropy.io import ascii
from astropy.time import Time, TimeDelta
from astropy.table import unique
import csv, requests, sys, json
import numpy as np

target_lists = {
    'nickel': 'https://ziggy.ucolick.org/yse/explorer/10/download?format=csv',
    'lcogt': 'https://ziggy.ucolick.org/yse/explorer/13/download?format=csv'
}

display = False

def download_targets(telescope):

    # Resolve the URL for the correct target list
    if telescope.lower() in target_lists.keys():
        url = target_lists[telescope.lower()]

        # Use requests to get list of targets
        data = requests.get(url)

        # Format into a table with the same names as standard target file
        table = ascii.read(data.text)
        table.sort('obs_date')
        table = unique(table, keys='name', keep='last')

        return(table)

    # Can't resolve list so throw print an error and return None
    else:
        error = 'ERROR: could not resolve a target list for telescope={tel}'
        print(error.format(tel=telescope.lower()))
        return(None)

# Instantiate a lcogt object with the correct shibboleth and download the list
# of pending and completed request groups obtained recently
lco = lcogt('/home/ckilpatrick/scripts/shibboleth')
sdate = Time(datetime.now()) - TimeDelta(28, format='jd')
rg = lco.get_requestgroups(propid=lco.proposals, itype='1M0-SCICAM-SINISTRO',
    sdate=sdate)

# Get all of the targets for which we might want LCOGT imaging
print('Grabbing targets from YSE PZ...')
targets = download_targets('lcogt')

if not targets or len(targets)==0:
    print('ERROR: could not get any LCOGT targets from YSE PZ.')
    sys.exit()

print('#'*80)
print('#'*80)
print('There are {n} targets'.format(n=len(targets)))
print(targets['name'].data)
print('\n\n')


if display:
    for g in rg:
        for r in g['requests']:
            for w in r['windows']:
                print(g['name'],
                    r['configurations'][0]['instrument_type'], w['start'],
                    g['state'])

for target in targets:
    needs_obs = True
    cadence = lco.params['strategy']['default']['cadence']
    targcoord = SkyCoord(target['ra'], target['dec'], unit='deg')
    mag = target['min_mag']+0.03*(Time(datetime.now()).mjd -\
        Time(target['min_date']).mjd)

    print(lco.params['strategy']['default']['ipp'])
    lco.params['strategy']['default']['ipp']=1.0
    lco.params['strategy']['default']['window']=1.0

    now = Time(datetime.now())+TimeDelta(7*3600*u.s)
    now = now.mjd

    #if '2021dov' in target['name'] and Time('2021-02-22').mjd-now<5:
    #    cadence = 0.4
    #    lco.params['strategy']['default']['window']=0.4
    #    lco.params['strategy']['default']['ipp']=1.06
    if '2021hiz' in target['name'] and np.abs(Time('2021-03-30').mjd-now)<10:
        lco.params['strategy']['default']['ipp']=1.12
        lco.params['strategy']['default']['window']=0.16
        cadence = 0.16

    if '2021lxb' in target['name'] and np.abs(Time('2021-05-15').mjd-now)<5:
        lco.params['strategy']['default']['ipp']=1.4
        lco.params['strategy']['default']['window']=0.06
        cadence = 0.06

    print('\n\n################################')
    print('Analyzing scheduling of:',target['name'])

    # HACK - TODO FIX THIS
    for request in rg:
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
        now = Time(datetime.now()).mjd
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

    # Message about target status
    if needs_obs:
        message = '{target} needs a new observation.  '
        message += 'Calculating observation properties...'
        print(message.format(target=target['name']))

        ra = targcoord.ra.degree
        dec = targcoord.dec.degree

        propidx=0
        try_next_proposal = True

        while try_next_proposal:

            response = lco.make_obs_request(target['name'], ra, dec,
                mag, propidx=propidx)

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
                print(message.format(target=target['name']))
            else:
                message = 'Successfully scheduled {target}'
                print(message.format(target=target['name']))
        else:
            message = '{target} exposures not possible with current settings'
            print(message.format(target=target['name']))
            print(response)
    else:
        message = '{target} does not need a new observation.'
        print(message.format(target=target['name']))

#!/usr/bin/env python3

import lcogt, json, sys, os, requests, time
import numpy as np, glob, tarfile, shutil
from astropy.time import Time, TimeDelta
from astropy.coordinates import SkyCoord
from astropy import units as u
from astropy.io import ascii, fits
from datetime import datetime
import subprocess
#import astrodash

# Options - these can/should be modified if the script is run by someone else
# or on a different machine
shibboleth = '/home/ckilpatrick/scripts/shibboleth'
floydscmd = 'floydsauto --site {site} --input-root-dir {dir} '+\
            ' --output-root-dir {dir} -A'
stagedir = '/data2/LCOGT/FLOYDS/'
proposals = ['NOAO2019A-020','NOAO2017AB-012','NOAO2018A-005',
                'NOAO2019A-020-TC','NOAO2017AB-005','NOAO2018A-007',
                'NOAO2018B-022','NOAO2019A-001','NOAO2018B-022b',
                'NOAO2019B-008','NOAO2019B-004','KEY2017AB-001',
                'NOAO2019B-009','NAOC2017AB-001','LCO2018A-004',
                'FTPEPO2014A-004','ARI2017AB-002','NOAO2020A-008',
                'NOAO2020A-012','NOAO2020B-009','NOAO2020B-011',
                'NOAO2021A-016','NOAO2021A-001']
datetime_fmt = '%Y-%m-%d %H:%M:%S'
trim_range = [3200.0, 10500.0]

def yse_upload(target):
    command = '/home/ckilpatrick/scripts/python/mkFakeEmails.py '
    command += '-s '+target

    os.system(command)

# Downloads an essentially complete list (look at SQL query params) of all
# supernovae in the YSE PZ data base
def download_yse_pz_supernovae():
    url = 'https://ziggy.ucolick.org/yse/explorer/45/download?format=csv'
    data = requests.get(url)
    table = ascii.read(data.text)

    for key in table.keys():
        table.rename_column(key, key.lower())

    for key in table.keys():
        if 'transient_ra' in key: table.rename(key, 'ra')
        if 'transient_dec' in key: table.rename(key, 'dec')

    return(table)

def upload_spectrum_to_yse(file, inst):
    script = '/home/marley/photpipe/pythonscripts/uploadTransientData.py'
    args = ['-i',file,'--clobber','--obsgroup','UCSC',
        '--permissionsgroup','UCSC','--instrument',inst,
        '-s','/home/ckilpatrick/scripts/settings.ini','-e','--spectrum']

    try:
        r = subprocess.run([script]+args, timeout=20)
    except subprocess.TimeoutExpired as e:
        print('YSE PZ upload for',file,'expired')

def clean_up(directory):
    clean_files = ['arc*', 'sens*', 'tt*', 'flat*', '*blue*.fits',
        '*red*.fits', 'logfile', 'README']

    for pattern in clean_files:
        files = glob.glob(directory + '/' + pattern)
        for file in files:
            os.remove(file)

    subdirs = [directory + '/database', directory + '/pyraf']
    for subdir in subdirs:
        if os.path.exists(subdir):
            shutil.rmtree(subdir)

def reformat_sn_name(obj):
    field = ''
    if obj.upper().startswith('20'):
        field = obj.lower()
    if (obj.upper().startswith('AT20') or obj.upper().startswith('SN20')
        or obj.upper().startswith('AT_20') or obj.upper().startswith('SN_20')
        or obj.upper().startswith('AT 20') or obj.upper().startswith('SN 20')):
        field = obj[2:].lower()
        field = field.replace('_', '')
        field = field.replace(' ', '')
    if len(field)==5:
        field = field.upper()

    return(field)


def format_spectrum(file, field=None, trim=True, run_dash=False, hostz=None):

    hdu = fits.open(file)
    path = os.path.split(file)[0]
    if not path: path = '.'

    header = hdu[0].header

    wvl = np.arange(header['XMIN'], header['XMAX'], header['CD1_1'])
    flux = hdu[0].data[0][0]

    if not field:
        field = reformat_sn_name(header['OBJECT'])
    inst = 'FLOYDS-N'
    if header['SITEID'] == 'coj': inst = 'FLOYDS-S'
    t = Time(header['DATE-OBS'])
    coord = SkyCoord(header['RA'], header['DEC'], unit=(u.hour, u.degree))
    newheader = {'RA': coord.ra.degree, 'DEC': coord.dec.degree,
        'INSTRUMENT': inst, 'OBS_DATE': t.datetime.strftime(datetime_fmt),
        'OBS_GROUP': 'UCSC', 'SNID': field, 'GROUPS': 'UCSC,YSE'}

    if run_dash:
        # Create a temporary file to hold the spectrum while we run astrodash
        tmpfile = path + '/tmp'
        with open(tmpfile, 'w') as t_file:
            for wave,flx in zip(wvl, flux):
                if trim:
                    if wave < trim_range[0] or wave > trim_range[1]:
                        continue
                t_file.write('%s %s \n' % (wave,flx))

        # Now run astrodash on tmpfile
        if hostz:
            classification = astrodash.Classify([tmpfile],[hostz],
                rlapScores=True, knownZ=True, smooth=6)
        else:
            classification = astrodash.Classify([tmpfile],rlapScores=True,
                smooth=6)
        bestFits,z,bestTypes,rlapFlags,matchesFlag,zerr = \
            classification.list_best_matches()

        bestType='' ; bestz = -99 ; rlap = -99 ; bestzerr = -99 ; zquality = 0
        data_quality = 0 ; phase = -99

        if len(bestTypes)>0:
            rlapDat = rlapFlags[0]
            if len(bestTypes)>0:
                bestType = bestTypes[0][0]
                # If we got a good fit, mark data quality as good
                if float(bestTypes[0][2]) > 0.5:
                    data_quality=1
                phase = np.mean([int(v) for v in bestTypes[0][1].split(' to ')])

            if hostz:
                if len(z)>0: bestz = z[0]
                bestzerr = 0.0
                zquality = 1
            else:
                if (z and len(z)>0): bestz = z[0]
                if (zerr and len(zerr)>0): bestzerr = zerr[0]
                if bestzerr > -99 and bestz/bestzerr < 0.2: zquality=1

            if 'rlap' in rlapDat: rlap=float(rlapDat.split('rlap:')[1])

        newheader['rlap'] = rlap
        newheader['redshift'] = bestz
        newheader['redshift_err'] = bestzerr
        newheader['redshift_quality'] = zquality
        newheader['spec_phase'] = phase
        newheader['spectrum_notes'] = matchesFlag[0]+', '+bestType
        newheader['data_quality'] = data_quality

        # Finally, delete tmpfile
        if os.path.exists(tmpfile):
            os.remove(tmpfile)

    newfile = field + '.' + t.datetime.strftime('%Y%m%d') + '.' + inst + '.flm'
    outfile = path + '/' + newfile
    with open(outfile, 'w') as s_file:
        s_file.write('# wavelength flux \n')
        for key in newheader.keys():
            line = '# {key} {val} \n'.format(key=key, val=newheader[key])
            s_file.write(line)

        for wave, flx in zip(wvl, flux):
            if trim:
                if wave < trim_range[0] or wave > trim_range[1]:
                    continue
            s_file.write('%s %s \n' % (wave,flx))

    return(outfile, newheader)

def run_floydsauto(directory):
    # Check for tar file, if it exists then skip floydsauto
    if len(glob.glob(directory + '/*/*/*/specproc/*.gz'))>0:
        print('WARNING: floydsauto already ran on dir={0}'.format(directory))
        print('Exiting...')
        return(1)
    else:
        if not os.path.exists(directory):
            os.makedirs(directory)
        os.chdir(directory)

        # Validate that current working dir is not home dir
        cwd = os.getcwd()
        if cwd == '/home/ckilpatrick':
            return(0)
        else:
            # First check which site and then format floydscmd
            if len(glob.glob('coj*'))>0:
                site = 'coj'
            elif len(glob.glob('ogg*'))>0:
                site = 'ogg'
            else:
                print('ERROR: could not determine site')
                print('Exiting...')
                return(0)

            cmd = floydscmd.format(site=site, dir=directory)
            print(cmd)
            os.system(cmd)

            return(0)

def unpack_spectrum(directory):
    # Check for ntt*merge*.fits file.  If it exists then skip
    os.chdir(directory)
    spectrum = directory + '/ntt*merge*.fits'
    if len(glob.glob(directory + '/ntt*merge*.fits'))>0:
        spectrum_file = glob.glob(directory + '/ntt*merge*.fits')[0]
        print('WARNING: spectrum file={0} already exists'.format(spectrum_file))
        print('Exiting...')
        return(spectrum_file)
    else:
        tarball = directory+'/*/*/*/specproc/*.gz'
        if len(glob.glob(tarball))>0:
            tar = tarfile.open(glob.glob(tarball)[0])
            tar.extractall()
            tar.close()

            spectrum_file = glob.glob(spectrum)
            if len(spectrum_file)==0:
                print('ERROR: floydsauto did not output a merged spectrum!')
                return(None)
            else:
                return(spectrum_file[0])
        else:
            message = 'ERROR: no tar.gz file in directory={dir}'
            print(message.format(dir=directory))
            return(None)

lco = lcogt.lcogt(shibboleth)
start = Time(sys.argv[1]) - TimeDelta(2, format='jd')
end = Time(sys.argv[1]) + TimeDelta(0.5, format='jd')
obslist = lco.get_obslist(sdate=start, edate=end,
    propid=lco.proposals, telid='2m0a', obstype='SPECTRUM', rlevel=0)

sn_list = download_yse_pz_supernovae()

download_list = []
print(obslist)
# Check obslist for objects that look like supernovae and try to find calibrations
for obs in obslist:

    field = reformat_sn_name(obs['OBJECT'])
    site = obs['SITEID']

    # Could be a supernova, check for match with YSE PZ
    if field:
        message = 'Checking YSE PZ for match with field={field}'
        print(message.format(field=field))

        c1 = SkyCoord(obs['area']['coordinates'][0][0][0],
            obs['area']['coordinates'][0][0][1], unit='deg')

        radius = 600. # arcsec
        sn_match = (sn_list['ra']-c1.ra.degree)**2 +\
            (sn_list['dec']-c1.dec.degree)**2 < (radius/3600)**2
        sn = sn_list[sn_match]
        print(sn)

        if len(sn)==0:
            continue

        if (field not in sn['name'].data and
            field.upper() not in sn['name'].data):
            print(field.upper(),'not in matched objects')
            continue

        sn_match = (sn['name']==field) | (sn['name']==field.upper())
        sn = sn[sn_match][0]

        print('Got match to {0}'.format(sn['name']))

        hostz = None
        if sn['redshift'] is not None:
            hostz = float(sn['redshift'])

        # Definitely want to download, first check for calibrations
        deltat = TimeDelta(0.1, format='jd')
        t = Time(obs['DATE_OBS'])

        # Get observations from the same request as the target (i.e.,
        # arc and flats)
        dllist = lco.get_obslist(sdate=t-deltat, edate=t+deltat,
                propid=lco.proposals, telid='2m0a',
                reqnum=obs['REQNUM'], rlevel=0)

        # Now we have download list, want to download.  Get outdir
        outdir = stagedir + t.datetime.strftime('ut%y%m%d') + '/' +\
                field + '/' + str(obs['REQNUM'])

        if not os.path.exists:
                os.makedirs(outdir)

        lco.download_obslist(dllist, outrootdir=outdir, use_basename=True,
                skip_header=True, funpack=False)

        # Check if there are ARC and LAMPFLAT frames in outdir
        check = False
        all_files = glob.glob(outdir+site+'2m*00.fits')
        obstypes = [fits.getval(file, 'OBSTYPE') for file in all_files]

        lco.get_spectral_calibrations(t, '2m0a', site,
                    outrootdir=outdir, funpack=False)

        dllist = lco.get_standardobs()
        #lco.download_obslist(dllist, outrootdir=outdir, use_basename=True,
        #    skip_header=True, funpack=False)

        run_floydsauto(outdir)
        spectrum_file = unpack_spectrum(outdir)
        if spectrum_file:
            formatted_spectrum, header = format_spectrum(spectrum_file,
                    run_dash=False, hostz=hostz)
            message = 'Final formatted file is {0}.  Uploading to YSE PZ.'
            print(message.format(formatted_spectrum))

            if formatted_spectrum:
                upload_spectrum_to_yse(formatted_spectrum, header['INSTRUMENT'])

        clean_up(outdir)

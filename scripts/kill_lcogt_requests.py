import requests
from urllib.parse import urlencode

PORTAL_URL='https://observe.lco.global'


def get_username_password():
    return('ckilpatrick@northwestern.edu','g6rfWJD@93QtaFH^bw$$')

def get_token_header(username, password):

    # Need to generate a new token
    data = {'username': username, 'password': password}
    uri = 'https://observe.lco.global/'+'api/api-token-auth/'

    # Now run a request
    response = requests.post(uri, data=data)
    response = response.json()

    # Check that the request worked
    if 'token' in response.keys():
        fmt =  'Token {token}'
        header = {'Authorization': fmt.format(token=response['token'])}
        return(header)

def get_requestgroups(headers, propid=['NSF2023A-011']):

    params = {}
    results = []
    propids = []
    if propid is not None:
        propids = list(propid)

    for pid in propids:
            params['proposal'] = pid
            params['limit'] = 500

            # Now do request
            response = requests.get('https://observe.lco.global/'+'api/requestgroups/',
                params=params, headers=headers).json()

            results += response['results']

    return(results)

def cancel_observation(observation_id, headers):
    requestgroup_id = _get_requestgroup_id(observation_id, headers)

    response = make_request(
            'POST',
            f'{PORTAL_URL}/api/requestgroups/{requestgroup_id}/cancel/',
            headers=headers
    )

    return response.json()['state'] == 'CANCELED'

def _get_requestgroup_id(observation_id, headers):
    query_params = urlencode({'request_id': observation_id})

    response = make_request(
            'GET',
            f'{PORTAL_URL}/api/requestgroups?{query_params}',
            headers=headers
    )
    requestgroups = response.json()

    if requestgroups['count'] == 1:
        return requestgroups['results'][0]['id']
    else:
        raise Exception('Could not get request group ID')

def make_request(*args, **kwargs):
    response = requests.request(*args, **kwargs)
    response.raise_for_status()
    return response

args = get_username_password()
headers = get_token_header(*args)
data = get_requestgroups(headers)
for d in data:
    for obs in d['requests']:
        if obs['state']!='PENDING': continue
        print(obs['id'])
        cancel_observation(obs['id'], headers)

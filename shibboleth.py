import os, sys
from astropy.io import ascii

class shibboleth(object):
    def __init__(self, shibboleth):

        self.shibboleth = None

        # Verify that only the owner can read this file
        st = os.stat(shibboleth)
        permissions = oct(st.st_mode)[-3:]
        if permissions!='440':
            message = 'ERROR: we only accept owner read-only permissions!'
            print(message)
            sys.exit(1)
        else:
            self.shibboleth = ascii.read(shibboleth,
                names=['type','username','password'])

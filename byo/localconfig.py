# -*- coding: utf-8 -*-
import os

system_root = "/data/BIO2/pcp/systems"
#system_root = "/home/mjens/pcp/systems_bak"
ucsc_server = "141.80.186.52"
www_server = "141.80.186.52"
www_root = "/data/BIO2/ucsc_static/internal_only"
url_root = "/cgi-bin/rbp_browser/"
rel_root = "/ucsc_static/internal_only"

#www_root = "/home/mjens/www/ucsc_static"
#www_server = "141.80.188.111"
#ucsc_server = "dorina.mdc-berlin.de"
#rel_root = "/ucsc_static"


systems = {
    'hg38' : os.path.join(system_root,'hg38'),
    'hg19' : os.path.join(system_root,'hg19'),
    'hg18' : os.path.join(system_root,'hg18'),
    'ce6' : os.path.join(system_root,'ce6'),
    'rn4' : os.path.join(system_root,'rn4'),
    'rn6' : os.path.join(system_root,'rn6'),
    'mm10' : os.path.join(system_root,'mm10'),
    'mm9' : os.path.join(system_root,'mm9'),
    'dm3' : os.path.join(system_root,'dm3'),
    #'smed1' : os.path.join(system_root,'smed1'),
    'sacCer3' : os.path.join(system_root,'sacCer3'),
    'danRer7' : os.path.join(system_root,'danRer7'),
}

default_system = 'hg19'


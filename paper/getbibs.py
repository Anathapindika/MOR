from urllib.request import urlopen
from sys import argv

bfil = open('many.bib','a')

def getbib(code):
    url = 'http://adsabs.harvard.edu/cgi-bin/nph-bib_query?bibcode=' + \
          code + '&data_type=BIBTEX'
    page = urlopen(url)
    s = 0
    for l in page:
        line = str(l,encoding='utf-8')
        if s < 5:
            s += 1
        else:
            bfil.write(line)
    print ('Read',code)

if len(argv) > 1:
    code = argv[1]
    getbib(code)


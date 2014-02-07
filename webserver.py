#!/usr/bin/env python

import sys
import BaseHTTPServer
from join_lobstr import read_lobstr
from urlparse import urlparse, parse_qs

def tups_to_table(tups, title='', style=''):
    cell_width = len(tups[0])
    s = '<table>\n'
    s += '<tr><th colspan="%s" style="font-weight: bold;%s">%s</th></tr>' % (cell_width,style,title)
    s += '\n'.join('<tr>' + ' '.join('<td>%s</td>' % str(x) for x in t) + '</tr>' for t in tups)
    #s += '\n'.join([ '<tr><td>%s</td><td>%s</tr>' % (k,v) for k,v in tups ])
    s += '</table>\n'
    return s
    

def read_sputnik(f):
    results = {}
    first_line = True
    for line in f:
        fields = map(str.strip, line.split('\t'))
        if first_line:
            colnames = fields
            first_line = False
            continue

        info = dict(zip(colnames, fields))
        results[info['chr'] + ':' + info['start']] = fields

    return (colnames, results)


def get_sputnik(chrom, start):
    """returns (HTML table describing the sputnik info, the end coordinate)"""
    cols, data = sputnik
    try:
        entry = data[str(chrom) + ':' + str(start)]
    except KeyError:
        return 'no entry in sputnik table for chr=%s, start=%s' % (chrom, start)

    d = dict(zip(cols, entry))
    return (tups_to_table(zip(cols, entry), 'SPUTNIK'), d['end'])


def get_lobstr_one(sample, chrom, start, end):
    """Return a list of (key, value) tuples or an empty list."""
    t = trees[sample]['chr'+chrom]
    hits = t.find(int(start), int(end))
    # Just show the first hit for now
    if hits:
        # stitched.depth is really long and has no whitespace to wrap
        # (tuples are immutable so we rebuild the tuple list)
        hits[0] = [ tup if tup[0] != 'stitched.depth' else (tup[0], tup[1].replace('/', ' ')) for tup in hits[0] ]

    return [('experiment', sample)] + hits[0] if hits else []


def get_lobstr(samples, chrom, start, end):
    lobstrs = [ get_lobstr_one(s, chrom, start, end) for s in samples ]
    lobstrs = filter(None, lobstrs)
    if not lobstrs:
        return 'No lobSTR matches.'

    # match is a list of (k,v) tuples; the k is the header for each tuple
    headers = [ [ k for k,v in match ] for match in lobstrs ][0]

    # merge all of the data into one big tuple
    data = [ [ v for k,v in match ] for match in lobstrs ]

    # Zip across the headers and data in each row
    table = zip(*([ headers ] + data))
    return tups_to_table(table, 'lobSTR')
 

class Handler(BaseHTTPServer.BaseHTTPRequestHandler):
    def do_GET(self):
        print('__________GOT REQUEST___________')
        params = parse_qs(urlparse(self.path).query)
        print(params)
        chrom, start = (-1, -1)
        html = 'unspecified error'
        try:
            chrom, start = (params['chr'][0], params['start'][0])
            html, end = get_sputnik(chrom, start)
            html += '<br/><hr/><br/>\n'
            html += get_lobstr(samples, chrom, start, end)
        except KeyError:
            html = 'ERROR: must pass chr and start arguments in the URL.<br/>'
            html += 'E.g., localhost:port/?chr=1&start=958201'

        self.__writeout(html)
        

    def __writeout(self, html_string):
        resp = "<html><body>%s</body></html>" % html_string
        self.send_response(200)
        self.send_header("Content-type", "text/html")
        self.send_header("Content-length", len(resp))
        self.end_headers()
        self.wfile.write(resp)
        self.wfile.close()

if len(sys.argv) != 2:
    print('usage: %s candidate_table' % sys.argv[0])
    exit(1)


samples = [ 'heart_bulk', 'cortex_bulk', 'neuron_2', 'neuron_3', 'neuron_51', 'neuron_6', 'neurons_100batch' ]
# for testing
#samples = [ 'heart_bulk' ] #, 'neuron_2' ]
print('loading SPUTNIK results..')
f = open(sys.argv[1], 'r')
sputnik = read_sputnik(f)

print('building lobSTR indexes for samples: ' + str(samples))
trees = read_lobstr(samples)


address = ('127.0.0.1', 3790)
print('starting webserver on ' + str(address))
httpd = BaseHTTPServer.HTTPServer(address, Handler)
while True:
    httpd.handle_request()

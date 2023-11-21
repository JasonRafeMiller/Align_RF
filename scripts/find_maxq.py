import sys
import gzip

SAMPLE_SIZE = 100000
NUM_LINES = 4 * SAMPLE_SIZE

filename = sys.argv[1]
print('Input file:', filename)

line_count = 0
lines1234 = 0
maxord = 0
maxchr = None
with gzip.open(filename,'r') as fin:
    for line in fin:
        lines1234 += 1
        line_count += 1
        if lines1234 == 4:
            lines1234 = 0
            this_maxord = max(line)
            if this_maxord > maxord:
                maxord = this_maxord
                maxchr = chr(maxord)
                print('Max so far:', maxchr, '=', maxord)
        if line_count > NUM_LINES:
            break
print('Num sequences examined:', SAMPLE_SIZE)
print('Num lines examined:', NUM_LINES)
print('The maximum QV was:', maxchr, '=', maxord)

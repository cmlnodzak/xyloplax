#!/usr/bin/py

# This script may be useful for searching a specific sequence against the blast 'nr' database. 
# In this case, the identifier is given by the record.id 'A_Xyloplax@1'. This can be modified for any ID in an MSA formatted file. 
# change the AlignIO.read 'format' argument given the filetype.

from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML
from Bio import AlignIO
import sys

filename = sys.argv[1]
outfile = str(filename)+'_blast.txt'

alignment = AlignIO.read(filename, 'phylip-relaxed')

with open(outfile, 'w') as blastout:
    for record in alignment:
        if record.id == 'A_Xyloplax@1':
            blast_record = NCBIWWW.qblast('blastp','nr', record.seq, matrix_name='BLOSUM62', gapcosts='11 1', word_size=6) ## NCBI default settings
            blast = NCBIXML.read(blast_record)
            for alignment in blast.alignments:
                for hsp in alignment.hsps:
                    if hsp.expect < 0.001: ### The 'alignment.title' species names are specific to find significant hits closely related to Asteroides
                        if 'Strongylocentrotus' in alignment.title:
                            blastout.write('TITLE: '+ str(alignment.title)+'\n')
                            blastout.write('LENGTH '+str(alignment.length)+'\n')
                            blastout.write(str(hsp.expect)+'\n'+'\n')
                        elif 'Saccoglossus' in alignment.title:
                            blastout.write('TITLE: '+ str(alignment.title)+'\n')
                            blastout.write('LENGTH: '+str(alignment.length)+'\n')
                            blastout.write('E-Val: '+str(hsp.expect)+'\n\n')
                        elif 'Patiria miniata' in alignment.title:
                            blastout.write('TITLE: '+ str(alignment.title)+'\n')
                            blastout.write('LENGTH: '+str(alignment.length)+'\n')
                            blastout.write('E-Val: '+str(hsp.expect)+'\n\n')
                        elif 'Asterias rubens' in alignment.title:
                            blastout.write('TITLE: '+ str(alignment.title)+'\n')
                            blastout.write('LENGTH: '+str(alignment.length)+'\n')
                            blastout.write('E-Val: '+str(hsp.expect)+'\n\n')
                        elif 'Heliocidaris crassi spina' in alignment.title:
                            blastout.write('TITLE: '+ str(alignment.title)+'\n')
                            blastout.write('LENGTH: '+str(alignment.length)+'\n')
                            blastout.write('E-Val: '+str(hsp.expect)+'\n\n')
                        elif 'Patiria pectinifera' in alignment.title:
                            blastout.write('TITLE: '+ str(alignment.title)+'\n')
                            blastout.write('LENGTH: '+str(alignment.length)+'\n')
                            blastout.write('E-Val: '+str(hsp.expect)+'\n\n')
                        elif 'Apostichopus japonicus' in alignment.title:
                            blastout.write('TITLE: '+ str(alignment.title)+'\n')
                            blastout.write('LENGTH: '+str(alignment.length)+'\n')
                            blastout.write('E-Val: '+str(hsp.expect)+'\n\n')
    print('BLAST search complete for: '+str(filename)+'\nPreparing next alignment file...')

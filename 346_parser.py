#!/usr/bin/py
from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML
from Bio import AlignIO
import sys
import operator

filename = sys.argv[1]
outfile = str(filename)+'_blast.txt'
print('Preparing alignment file: '+str(filename)+'...')
ids_list = []
lengths = []
stats = {}
alignment = AlignIO.read(filename, 'phylip-relaxed')
amino_chars = ['G', 'A', 'L', 'M', 'F', 'W', 'K', 'Q', 'E', 'S', 'P', 'V', 'I', 'C', 'Y', 'H', 'R', 'N', 'D', 'T']
with open(outfile, 'w') as blastout:
    for record in alignment:
        if record.id == 'Saccoglossus@1':
            continue
        elif record.id == 'refseqspurp@1':
            continue
        else:
            ids_list.append(record.id)
            seqlength = 0
            for letter in record.seq:
                if letter in amino_chars:
                    seqlength += 1
        lengths.append(seqlength)
    stats = dict(zip(ids_list, lengths))
    blast_seq = max(stats.iteritems(),key = operator.itemgetter(1))[0]
    blastout.write('The sequence blasted was: '+str(blast_seq)+'\n\n')
    print('The record blasted was: '+str(blast_seq))
    for record in alignment:
        if blast_seq == record.id:
            blast_record = NCBIWWW.qblast('blastp','nr', record.seq, matrix_name='BLOSUM62', gapcosts='11 1', word_size=6)
            blast = NCBIXML.read(blast_record)
            for alignment in blast.alignments:
                for hsp in alignment.hsps:
                    if hsp.expect < 0.001:
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


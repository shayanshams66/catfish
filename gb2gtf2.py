#!/usr/bin/env python
import os, sys
from datetime import datetime
from Bio      import SeqIO

def gb2gtf( source='gb2gtf',allowedTypes=set(['gene','CDS','tRNA','tmRNA','rRNA','ncRNA','mRNA']) ):
  """
  """
  handle = sys.stdin
  for gb in SeqIO.parse( handle,'gb' ):
    acc     = gb.name #gb.id #gb.name #gb.description # # 
    skipped = 0
    skippedTypes = set()
    for f in gb.features:
      #process only gene and CDS entries
      if f.type not in allowedTypes:
        skipped += 1
        skippedTypes.add( f.type )
        continue
      exom=''  
      comments = ''
      if 'locus_tag' in f.qualifiers:
        #use locul tag as gene_id/transcript_id
        gene_id = transcript_id = f.qualifiers['locus_tag'][0]
        comments = 'gene_id "%s"; transcript_id "%s"' % ( gene_id,transcript_id )
        
      elif 'gene' in f.qualifiers:
        gene_id = transcript_id = f.qualifiers['gene'][0]
        comments = 'gene_id "%s"; transcript_id "%s"' % ( gene_id,transcript_id )        

      if not comments:
        sys.stderr.write( "Error: Neither `gene` nor `locus_tag` found for entry: %s\n" % '; '.join( str(f).split('\n') ) )
        continue
              
      if   'product' in f.qualifiers:
        comments += '; protein_id "%s"' % f.qualifiers['product'][0]
      elif 'protein_id' in f.qualifiers:
        comments += '; protein_id "%s"' % f.qualifiers['protein_id'][0]

      #add external IDs  
      if 'db_xref' in f.qualifiers:
        for extData in f.qualifiers['db_xref']:
          comments += '; db_xref "%s"' % extData
      
      if int(f.strand)>0: strand = '+'
      else:               strand = '-'
      gtf = '%s\t%s\t%s\t%s\t%s\t%s\t.\t%s\t.\t%s' % ( acc,source,f.type,exom,f.location.start.position+1,f.location.end.position,strand,comments ) #f.frame,
      print gtf
      if len(f.location.parts)>1:
        for part in f.location.parts:
          comments = ''
          if 'locus_tag' in f.qualifiers:
            #use locul tag as gene_id/transcript_id
            gene_id = transcript_id = f.qualifiers['locus_tag'][0]
            comments = 'gene_id "%s"; transcript_id "%s"' % ( gene_id,transcript_id )
            
          elif 'gene' in f.qualifiers:
            gene_id = transcript_id = f.qualifiers['gene'][0]
            comments = 'gene_id "%s"; transcript_id "%s"' % ( gene_id,transcript_id )        

          if not comments:
            sys.stderr.write( "Error: Neither `gene` nor `locus_tag` found for entry: %s\n" % '; '.join( str(f).split('\n') ) )
            continue
                  
          if   'product' in f.qualifiers:
            comments += '; protein_id "%s"' % f.qualifiers['product'][0]
          elif 'protein_id' in f.qualifiers:
            comments += '; protein_id "%s"' % f.qualifiers['protein_id'][0]

          #add external IDs  
          if 'db_xref' in f.qualifiers:
            for extData in f.qualifiers['db_xref']:
              comments += '; db_xref "%s"' % extData
          
          #code strand as +/- (in genbank 1 or -1)
          if int(part.strand)>0: strand = '+'
          else:               strand = '-'
          gtf = '%s\t%s\texom\t%s\t%s\t.\t%s\t.\t%s' % ( acc,source,part.start.position+1,part.end.position,strand,comments ) #f.frame,
          print gtf
    sys.stderr.write( "%s\tSkipped %s entries having types: %s.\n" % ( gb.id,skipped,', '.join(skippedTypes) ) )

if __name__=='__main__': 
  t0=datetime.now()
  gb2gtf()
  dt=datetime.now()-t0
  sys.stderr.write( "#Time elapsed: %s\n" % dt )

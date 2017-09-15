# Commands from Sasson and Ryan, 2017


#### Align 18s sequences w SSU-Align version 0.1.1 (Nawrocki 2000):
```
export SSUALIGNDIR=/usr/local/share/ssu-align-0.1.1
ssu-align 18s.fa 18s.ssu
ssu-mask 18s.ssu
ssu-mask --stk2afa 18s.ssu
```
  
#### Trim the resulting alignment with Gblocks version 0.91b (Castresana 2000) using Gblockswrapper 
#### https://bitbucket.org/caseywdunn/labcode/src/master/scripts_phylogenomics_21Feb2009/Gblockswrapper
```
grep -v '^#' 18s.ssu/18s.ssu.eukarya.mask.stk | perl -ne 's/^/>/; s/\s+/\n/; print;' > 18s.aln.fa
perl -pi -e 's/^>\/\/\s*$//;' 18s.aln.fa
perl -pi -e 's/^>\s*$//;' 18s.aln.fa
Gblockswrapper 18s.aln.fa
```

#### Convert Gblocks output to FASTA format:
```
perl -pi -e 's/ //g' 18s.aln.fa-gb
fasta2phy.pl 18s.aln.fa-gb > 18s.aln.fa-gb.phy
```

#### Calculate 18s tree w RaxML version 8.1.21 constraining w composite tree “known_relationships.tre”:
```
raxmlHPC-PTHREADS-SSE3 -s 18s.aln.fa-gb.phy -n 18s.knownrel -m GTRCAT -T 20 -p 1234 -g composite.tre > rax.out 2> rax.err & 
```

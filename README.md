# Metaplot-and-Counting-fragments
Counting fragments and making metaplots using bedtools and python.

Author: Mrutyunjaya Parida, David Price Lab, UIOWA

## Counting fragments for given genomic intervals:

```bedtools coverage –a genomic-intervals.bed –b mapped.deduped.fragments.bed –counts > genomic-intervals-fragment-counts.bed```

## Counting sense fragments for given genomic intervals:

```bedtools coverage –a genomic-intervals.bed –b mapped.deduped.fragments.bed –s –counts > genomic-intervals-fragment-counts-sense.bed```
The –s option requires same strandedness.

## Counting anti-sense fragments for given genomic intervals:

```bedtools coverage –a genomic-intervals.bed –b mapped.deduped.fragments.bed –S –counts > genomic-intervals-fragment-counts-antisense.bed```
The –S option requires same strandedness.

## Counting 5’ ends of sense fragments for given genomic intervals:

First run this python code to get the 5’ ends of fragments: 
```
OFILE=open('mapped.deduped.fragments.5prime.bed', 'w'); #change name of the output file here
for line in open('mapped.deduped.fragments.bed','r'): #change name of the input file here
      DATA=line.strip().split('\t');
      if DATA[5]== '+':
          OFILE.write('\t'.join([DATA[0]]+[repr(int(DATA[1])), repr(int(DATA[1])+1)]+DATA[3:]).strip()+'\n');
      elif DATA[5]== '-':
          OFILE.write('\t'.join([DATA[0]]+[repr(int(DATA[2])-1),repr(int(DATA[2]))]+DATA[3:]).strip()+'\n');
OFILE.close();
```
Next run the following command:

```bedtools coverage –a genomic-intervals.bed –b mapped.deduped.fragments.5prime.bed –s –counts > genomic-intervals-fragment-5primeends-counts-sense.bed```

## Counting 5’ ends of anti-sense fragments for given genomic intervals:

```bedtools coverage –a genomic-intervals.bed –b mapped.deduped.fragments.5prime.bed –S –counts > genomic-intervals-fragment-5primeends-counts-anti-sense.bed```

## Counting 3’ ends of sense fragments for given genomic intervals:

First run this python code to get the 3’ ends of fragments: 
```
OFILE=open('mapped.deduped.fragments.3prime.bed', 'w'); #change name of the output file here
for line in open('mapped.deduped.fragments.bed','r'): #change name of the input file here
      DATA=line.strip().split('\t');
      if DATA[5]== '-':
          OFILE.write('\t'.join([DATA[0]]+[repr(int(DATA[1])), repr(int(DATA[1])+1)]+DATA[3:]).strip()+'\n');
      elif DATA[5]== '+':
          OFILE.write('\t'.join([DATA[0]]+[repr(int(DATA[2])-1),repr(int(DATA[2]))]+DATA[3:]).strip()+'\n');
OFILE.close();
```
Next run the following command:

```bedtools coverage –a genomic-intervals.bed –b mapped.deduped.fragments.3prime.bed –s –counts > genomic-intervals-fragment-3primeends-counts-sense.bed```

## Counting 3’ ends of anti-sense fragments for given genomic intervals:

```bedtools coverage –a genomic-intervals.bed –b mapped.deduped.fragments.3prime.bed –S –counts > genomic-intervals-fragment-3primeends-counts-anti-sense.bed```

## Counting depth for every position of given genomic intervals:

First, make sure all the genomic intervals are of the same length.
Replace the –counts option with the –d option and run the commands above for specific tasks. Do not use –s and –S if fragments or 5’ ends or 3’ ends from both strands are counted for a list of genomic intervals.
Metaplot after counting depth for every position of given genomic intervals.
First, make sure all the genomic intervals are of the same length.
Run this python code to get the summary (mean or sum) of fragment counts of all positions in the list of genomic intervals:
```
P={};N={};GenIntvls=set();Meta={};

def metaplotmean(Pos,Neg,P,N,Dirctn):
    Meta={};tmp=[];
    #change name of the output file here
    OFILE=open('mapped.fragments-mean.metaplot','w');
    for p in Pos:
        Meta.setdefault(p,0.0);
        Meta[p]+=float(P[p])/GenIntvls;

    for n in Neg:
        tmp.append(float(N[n])/GenIntvls);

    for pp in Pos:
        Meta[pp]+=tmp[pp-1];

    for o in Pos:
        if Dirctn=='sense':
            OFILE.write(repr(o).strip()+'\t'+repr(round(Meta[o],3)).strip()+'\n');
        elif Dirctn=='antisense':
            OFILE.write(repr(o).strip()+'\t'+repr(-1*round(Meta[o],3)).strip()+'\n');
    OFILE.close();

def metaplotsum(Pos,Neg,P,N,Dirctn):
    Meta={};tmp=[];
    #change name of the output file here
    OFILE=open('mapped.fragments-sum.metaplot','w');
    for p in Pos:
        Meta.setdefault(p,0);
        Meta[p]+=int(P[p]);

    for n in Neg:
        tmp.append(int(N[n]));

    for pp in Pos:
        Meta[pp]+=tmp[pp-1];

    for o in Pos:
        if Dirctn=='sense':
            OFILE.write(repr(o).strip()+'\t'+repr(round(Meta[o],3)).strip()+'\n');
        elif Dirctn=='antisense':
            OFILE.write(repr(o).strip()+'\t'+repr(-1*round(Meta[o],3)).strip()+'\n');
    OFILE.close();
    
#change name of the input file here
for line in open('mapped.fragments.depth','r'):
    DATA=line.strip().split('\t');
    GenIntvls.add('\t'.join(DATA[:6]).strip());
    if DATA[5].strip()=="+":
       P.setdefault(int(DATA[6]),0);
       P[int(DATA[6])]=P[int(DATA[6])]+int(DATA[7]);
    elif DATA[5].strip()=="-":
       N.setdefault(int(DATA[6]),0);
       N[int(DATA[6])]=N[int(DATA[6])]+int(DATA[7]);

GenIntvls=len(GenIntvls);
Pos=P.keys();
Pos.sort(key=int);

Neg=N.keys();
Neg.sort(key=int,reverse=True);

#When computing average for the sense data
metaplotmean(Pos,Neg,P,N,'sense');

#When computing average for the antisense data
metaplotmean(Pos,Neg,P,N,'antisense');

#When computing sum for the sense data
metaplotsum(Pos,Neg,P,N,'sense');

#When computing average for the antisense data
metaplotsum(Pos,Neg,P,N,'antisense');
```
## Output:
The metaplot data is present in the mapped.fragments-mean.metaplot and/or mapped.fragments-sum.metaplot files in a 2 column format. The first column shows hypothetical genomic positions of all genomic intervals used in these metaplots. The second column shows the summary data for all genomic positions. The metaplot data can be copied to an MS Excel file for plotting purposes.

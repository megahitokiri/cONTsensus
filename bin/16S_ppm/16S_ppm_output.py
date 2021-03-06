import re, os,sys

def filter_multiID(nome):
    g=open(nome).readlines()
    d1={}
    l1=[]
    for j in range(len(g)):
        i=g[j]
        if len(i) > 0:
            i2= i[:-1]
            Id=i2.split()[0]
            l1.append(Id)
            annota_l=i2.split("\t")[1:]
            annot=""
            for e in annota_l:
                annot=annot+"\t"+e
            if Id in d1:
               pass  
            else:
               d1[Id]= annot
    d2={}
    out=open(nome+"_best","w")
    for i in d1:

        if l1.count(i)> 1:
          d2[i]=1  

        out.write(i+d1[i]+"\n")
    out.close()


def filter_blast(top_hits, k):
    hits=open(top_hits).read().split('\n')
    filt_out=''
    for e in hits:
      if len(e) > 1:
       if len(e.split('\t')):
        cov=int(e.split('\t')[2])
        if cov > int(k):
            filt_out=filt_out+e+'\n'
    return filt_out


def filter_blast_percid(top_hits, p):
    hits=open(top_hits).read().split('\n')
    filt_out=''
    for e in hits:
      if len(e) > 1:
       if len(e.split('\t')):
        perc=float(e.split('\t')[-1])
        if perc >= float(p):
            filt_out=filt_out+e+'\n'
    return filt_out


def count_species(name):
    outname=name.replace("strain","species")
    out=open(outname,"w")
    h=open(name).read().split("\n")
    diz_sp={}
    for i0 in h:
       if len(i0)>0:
        n=i0.split()[0]
        if "Candidatus" in i0.split()[1]:
           i0=i0.replace("Candidatus ","Candidatus-")        
        k=i0.split()[1].replace(" ","")
        k2=i0.split()[1]+"_"+i0.split()[2]
        i=k2
        if i in diz_sp:
            diz_sp[i]=diz_sp[i]+int(n)
        else:
            diz_sp[i]=int(n)
    for i in diz_sp:
        out.write(i+"\t"+str(diz_sp[i])+"\n")

    out.close()

evalue=0.00001
align_coverage_cutoff=0
align_perc_id_cutoff=0

###evalue: evalue cut of for blastn
###align_coverage_cutoff: filter alignments for coverage above cutoff (e.g. if =60, filters out align <60% coverage)
###align_perc_id_cutoff: filter alignments for %identity > cutoff (e.g. if =70, filters out align <70% identity)

inpQ=sys.argv[1] 
inp2=inpQ.split(".unfilter_top_hits.txt")[0];inp2=inp2.split("/")[-1]

input_folder=sys.argv[2]
out_folder = sys.argv[3]

nome_temp=input_folder+inp2+".unfilter_top_hits.txt"
f2=filter_multiID(nome_temp)

koff=align_coverage_cutoff
p=align_perc_id_cutoff

print('Filtering best blast top hit candidates on:',nome_temp)
f=filter_blast(nome_temp+"_best",koff)
outo=open("cov_"+str(koff)+"_"+inp2+"_top_hits.txt","w")
outo.write(f)
outo.close()
print('Filtering by Percentage ID on top hit candidates on',nome_temp)
f=filter_blast_percid("cov_"+str(koff)+"_"+inp2+"_top_hits.txt",p)
outo=open("cov_"+str(koff)+"_idfilt_"+str(p)+"_"+inp2+"_top_hits.txt","w")
outo.write(f)
outo.close()
os.remove("cov_"+str(koff)+"_"+inp2+"_top_hits.txt")

os.system("mv cov_"+str(koff)+"_idfilt_"+str(p)+"_"+inp2+"_top_hits.txt "+nome_temp)

log=open("log_run_fromfasta_pipe_16S_NCBI.log","w")
log.write("Finished Filtering best hits in blast")

if koff > 0 or p> 0:
        log.write("\nBlast filtered for Coverage="+str(koff)+ " %identity="+str(p)+"\n")

print('Starting Strains Count on blast filtered file')		
os.system("cat "+str(nome_temp)+" | cut -f4 | sort | uniq -c | sort -nr > "+str(inp2)+"_strain_counts.txt")
nome_x_demo=str(inp2)+"_strain_counts.txt"
count_species(nome_x_demo)

ss=nome_x_demo.replace("strain","species")

os.system("sort -k 2 -nr "+ss+" > "+ss+".sorted")


inpo=open(ss+".sorted").readlines()

outo_name=("Out16S_"+ss)
outo=open(outo_name,"w")
num_reads=0
for m in inpo:
            m=m[:-1]
            m2=m.split()
            nm=int(m2[-1])
            num_reads=num_reads+nm
sample_name=inp2

outo.write("#Identification: "+sample_name+", reads="+str(num_reads)+"\n#Genus,species,counts,relative abundance (%)\n")
print('Writing output file')
for m in inpo:
            m=m[:-1]
            num=int(m.split()[-1])
        
            rel=num/float(num_reads)
            rel1=round(rel,8)
            rel1=100*rel1
            m2=m.split()
            line=m2[0]+","+str(m2[-1])+","+str(rel1)
            line=line.replace("_",",")
            outo.write(line+"\n")

outo.close()

### below, clean up of temporary files, and move output files to output folder

out=out_folder

if os.path.isdir(out):
    os.system("ls | wc -l")
    pass
else:
    os.makedirs(out)

log.close()
print('Cleaning Temporary files')
os.system("mv "+outo_name+" "+out)

os.system("mv log_run_fromfasta_pipe_16S_NCBI.log "+out)
os.system("rm "+ss+" "+ss+".sorted "+nome_x_demo+" "+nome_temp)

os.system("mv "+nome_temp+"_best"+" "+out)
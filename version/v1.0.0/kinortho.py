#-----------------------------------------------------------------------------
#----------------------------NEEDS TO RUN WITH LINUX--------------------------
#-----------------------------------------------------------------------------
import sys
import Bio          #------------ importing stuff----------------------------
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
import numpy as np
import os
import math
import statistics
#--------------------------------NEED BIOPYTHON TO RUN!-----------------------
#--------------------------------NEED BLAST+ TO RUN !-------------------------
#--------------------------------NEED NUMPY to RUN !!!!!!_--------------------


#-----------------------------------Helper varibal/input stuff----------------
infolder = sys.argv[1]
x=2
fullseq=""
domseq=""
out=""
evalue=""
thread=""
inflat=""
minev=""

#woop="1234567"
#woop=woop[0:8]
#print (woop)

#-----------------just taking the input from bash and making------------------
#-----------------it be able to take any comand in any order------------------

for i in range(int(len(sys.argv)/2)-1):#loops in pairs through the bash command
    if len(sys.argv)>x:                #       if the commands are even present
        if sys.argv[x]=="-f":
            fullseq = sys.argv[x+1]
        if sys.argv[x]=="-d":
            domseq = sys.argv[x+1]
        if sys.argv[x]=="-o":
            out = sys.argv[x+1]
        if sys.argv[x]=="-E":
            evalue = sys.argv[x+1]  #         all these conditionals check what
        if sys.argv[x]=="-t":       #              are present in the bash line
            thread = sys.argv[x+1]
        if sys.argv[x]=="-i":
            inflat = sys.argv[x+1]
        if sys.argv[x]=="-e":
            minev = sys.argv[x+1]
    x=x+2
if minev == "":
    minev = float(1e-200)
if evalue == "":
    evalue=float(1e-5)
if inflat== "":
    inflat=1.5
#firstcommand= "makeblastdb -in ref.fasta -dbtype prot -out ref"

#-------------------------------defining some methods/funtions-----------------
def build_blast_db(input,dbname):
    firstcommand= "makeblastdb -in "
    firstcommand= firstcommand + input
    firstcommand= firstcommand +" -dbtype prot -out "
    firstcommand= firstcommand + dbname
    os.system(firstcommand)


def homology_search(inf,db,outf):
    command="blastp -outfmt \"6 qseqid sseqid sstart send evalue bitscore positive\""
    command= command + " -query " + inf
    command= command + " -db " + db
    if(thread!=""):
        command= command + " -num_threads " + thread

    if(evalue!=""):
        command= command + " -evalue " + str(evalue)

    # add condtionals
    command= command + " > " + outf
    os.system(command)


class Edgy:
    def __init__(self,location,weight):
        self.w = weight
        self.loc = location

    def getW(self):
        return self.w

    def getL(self):
        return self.loc

#--------------------------------------------------------------------------
fullrec=[]
index=0  # helper varibal
rec=[]
tax=[]
fquary=[]
amountofspec=0
if fullseq != "":
    records = list(SeqIO.parse(fullseq, "fasta"))

    for record in records:
        r=record.id

        #r=r+'#' + infile
        record.id = r
        fullrec.append(r)
        index= index +1
        quary=""
        b1=0
        b2=0
        for char in r:
            if char == '|' and b1==1:
                b2=1
            if b1==1 and b2 ==0:
                quary=quary+char


            if char=='|':
                b1=1
        #print(record.id)
        fquary.append(quary)


#print(fullrec)
dquary=[]
if domseq !="":
    records = list(SeqIO.parse(domseq, "fasta"))
    for record in records:
        r=record.id
        quary=""
        b1=0
        b2=0
        for char in r:
            if char=='.':
                b1=1
            if b1==0:
                quary=quary+char


        #print(quary)
        dquary.append(quary)





#-------------------------This is for retreaving files------------------------
for infile in os.listdir(infolder):
    amountofspec=amountofspec+1
    goodfile = infolder + infile
    records = list(SeqIO.parse(goodfile, "fasta"))      # reads the input files
    for record in records:
        r=record.id
        r=r+'#' + infile
        record.id = r
        rec.append(record)
        index= index +1
print(index)
#---------------------------This is for writing a file------------------------
ref=open("ref.fasta","w+")
SeqIO.write(rec,"ref.fasta","fasta")
ref.close()







#woop=Edgy("asdklfal",3)
#print(woop.getW())
#print(woop.getL())


#'''

#-------------First step in the program to build a database--------------------
build_blast_db("ref.fasta","firstDB")
#---------------Blastp the full seq agasit the database------------------------
#'''
#'''

if fullseq != "":
    homology_search(fullseq,"firstDB","txtHelper.txt")

    print("first homology search done")


    helperfile=open("txtHelper.txt","r")
    lines = helperfile.readlines()

    #-----------------------read the results form the Blastp-----------------------
    helperfile=open("txtHelper.txt","r")
    lines = helperfile.readlines()

    sseqids=[]
    tax={}
    fulltax={}




    curfull="woop"
    # reading recording speiceis for each protien assighning a dictionary
    for line in lines:
        b=0
        b2=0
        fcur=""
        fkey=""
        ftax=""
        for char in line:
            if char=='\t':
                fkey=fcur
            if b == 0:
                if char != '\t':
                    fcur=fcur+char
            else:
                if char=='#':
                    b2=1

                if b2==1:
                    if char=='\t':
                        break
                    if char != '\t':
                        ftax=ftax+char
            if char=='\t':
                b=1
        if curfull!=fkey:
            curfull=fkey
            fulltax[fkey]=ftax


    #print(fulltax)

    index=0

    for line in lines:  #                for each line in the result in the blastp
        l=list(line)
        index=0
        tempstring=""
        helperbool=0
        b=0
        for chara in l:  #                    for each chararicter in the line
            if chara == '\t':
                for i in range(len(l)):#        record teh line until the tab char
                    if l[index+i+1] != '\t':
                        tempstring= tempstring + l[index+i+1]
                    else:
                        for sseqid in sseqids: #                   for each sseqid
                            if str(sseqid) == tempstring:
                                helperbool=1
                        if helperbool ==0:
                            sseqids.append(tempstring)#add ids to the array of ids
                            #print(tempstring)
                        break
                break
            index=index+1



    #-----------------------geting the files ready for all vs all blast------------

    newrec=[]
    for record in rec:


        for sseqid in sseqids:
            if record.id==sseqid:
                newrec.append(record)
    #print(newrec)

    print(len(rec))
    print(len(sseqids))
    print(len(newrec))

    #----------------read and writes the file for the second blastP---------------
    newref=open("homology.fasta","w+")
    SeqIO.write(newrec,"homology.fasta","fasta")
    newref.close()

    build_blast_db("homology.fasta","homology")

    #----------------blastp of all vs all----------------------------------------
    print("allvsall start")

    homology_search("homology.fasta","homology","all_vs_all.blastp")


    print("allvsall finish")


    #----------------------------------------------------------------------------
    #----------reads and retrives the data from the output of the----------------
    #------------all vs all blastP and puts into dictionaries--------------------
    #----------------------------------------------------------------------------

    helperfile=open("all_vs_all.blastp","r")
    lines = helperfile.readlines()

    keys=[]
    halfkeys=[]
    indhalfkeys=[]
    fkeysE={}
    keysB={}
    keysP={}
    minE={}
    maxB={}
    maxP={}
    tax={}
    seckeys={}
    firkeys={}
    samebool=[]
    taxorder={}
    last="woop"

    #sstart={}
    #send={}
    taxorderhelper=0

    #aveP={}
    #-------------go through the files and find all the info of the -------------
    #--------------------AllvsAll Blastp outputfile------------------------------

    #          things you should know about about blastp files
    #          there differnt info is sesperated by a tab (\t)
    #       in the code when you see tab(\t) its the end thing inbetween info


    for line in lines:# for each line
        l=list(line)
        #print(l)
        #helper varibles
        index=0
        fin =0
        f=""
        s=""
        firststring=""
        tempstring=""
        helperbool=0
        taxbool=0
        taxa=""

        for chara in l: #                      for each character in each line
            if chara !='\t':#     if not the end of the line add to temp sting
                tempstring= tempstring +chara
                if helperbool == 0:#          find first string and add to var
                    firststring= firststring+ chara
                    f=f+chara
                else:
                    s=s+chara
                #print(tempstring)
            if chara == '\t':#   if end of string in file add taxonomy to dic
                taxbool=0
                tax[tempstring]=taxa
                if not taxa in taxorder.keys():#       checks for no reapeats
                    taxorder[taxa]=taxorderhelper
                    taxorderhelper=taxorderhelper+1
                if helperbool==1:#               if the second id was passed
                    if f == s: #      if first id and sencond id is the same
                        samebool.append(1)#to be honest i dont even use this
                    else: #         but to scared to scared to get rid of it
                        samebool.append(0)
                    index=index+1       #index helps counts the where we are
                    keys.append(tempstring) #    add the full key to the dic

                    charbool=0              #       We are done with id part
                    secstring=""
                    for char in tempstring: #   for the charartors in stirng
                        if charbool==1: #    if on second id place it in dic
                            secstring=secstring+char
                        if char=='\t':
                            charbool=1
                    firkeys[tempstring]=firststring
                    seckeys[tempstring]=secstring

                    if not firststring in indhalfkeys:#checking for repeats
                        indhalfkeys.append(firststring)
    #                      make a array of half keys to help later funtions
    #                                       Both unquie and non unquie list
                    halfkeys.append(firststring)
                    #if not tempstring in aveB.keys():
                    #    aveB[tempstring]=[]
                    #    aveP[tempstring]=[]

            #------ after the ID and into the other data------
                    #ss=""
                    #if secstring not in sstart.keys():
                    #    sstart[secstring]=[]
                    for i in range(len(l)):# the start info and redoring it
                        if l[index]!='\t':
                    #        ss=ss+l[index]
                            index=index+1
                        if l[index]=='\t':
                            index=index+1
                            break
                    #if ss not in sstart[secstring]:
                    #sstart[secstring].append(ss)
                    #se=""
                    #if secstring not in send.keys():
                    #    send[secstring]=[]
                    for i in range(len(l)):#  the end info and recording it
                        if l[index]!='\t':
                    #        se=se+l[index]
                            index=index+1
                        if l[index]=='\t':
                            index=index+1
                            break
                    #if se not in send[secstring]:
                    #send[secstring].append(se)
                    e=""#                      the E-value and recording it
                    for i in range(len(l)):
                        if l[index]!='\t':
                            e=e+l[index]
                            index=index+1
                        if l[index]=='\t':
                            fkeysE[tempstring]=e
                            if not firststring in minE.keys():
                                minE[firststring]=e
                            else:
                                if float(e) < float(minE[firststring]):
                                    minE[firststring]=e
                            index=index+1
                            break
                    b=""#                   the bit-score and recording it
                    for i in range(len(l)):
                        if l[index]!='\t':
                            b=b+l[index]
                            index=index+1
                        if l[index]=='\t':
                            keysB[tempstring]=b
                            if not firststring in maxB.keys():
                                maxB[firststring]=b
                            else:
                                if float(b) > float(maxB[firststring]):
                                    maxB[firststring]=b

                            index=index+1
                            break
                    p=""#                    The postive and recording it
                    for i in range(len(l)):#   for the length of the line
                        if l[index]!='\t':#            tab inbetween info
                            p=p+l[index]
                            index=index+1 #             add more to index
                        if index==len(l):#                 if end of line
                            keysP[tempstring]=p#
                            if not firststring in maxP.keys():
                                maxP[firststring]=p#       add pos to dic
                            else:
                                if float(p) > float(maxP[firststring]):
                                    maxP[firststring]=p
                            fin=1
                            # helper bool that tells us that we are done
                            index=index+1 #          gotta add the index
                            break
                tempstring=tempstring+'\t'
                helperbool=1
                taxa=taxa+'#'

            if fin == 1:#                                  if done break
                break

            if taxbool==1:     #                   stuff for the tax dic
                taxa=taxa+chara
            if chara=='#':
                taxbool =1
            index=index+1

    helperfile.close()


    #print(len(keys))
    #for key in keys:
    #    print(key)
    firkeysf=firkeys
    seckeysf=seckeys

    print("allvall scanner all done")
    #----------------------creats helpers comp sturctures-----------------------

    tax2={}#                               helps dics for the taxonomy
    tax1={}
    taxorder={}
    taxorderhelper=0

    for key in keys: #                                for the all keys
        t1=""#                                               first tax
        t2=""#                                              second tax
        taxa=tax[key]#      gets the tax dic and slits it into two dic
        helperbool=0
        for char in taxa:
            if char=='#':
                helperbool=1
                continue
            if helperbool == 0:
                t1=t1+char
            if helperbool ==1:
                t2=t2+char
        if not t2 in taxorder.keys():
            taxorder[t2]=taxorderhelper
            taxorderhelper=taxorderhelper+1
        tax1[key]=t1
        tax2[key]=t2

    ftax1=tax1
    ftax2=tax2

    #-----------------------finding paralogs------------------------------------
    para={}
    p=[]
    i=0
    cur="woop"
    helperbool=0
    helperbool2=0
    for key in keys:

        if cur != halfkeys[i]:
            cur=halfkeys[i]
            helperbool=0
            para[halfkeys[i]]=[]
            helperbool2=0


        if tax2[key]!=tax1[key]:
            helperbool=1
        if helperbool == 0:
            if helperbool2==1:
                #print(key)
                #print(keysE[key])
                if float(fkeysE[key]) <= evalue:
                    if key not in p:
                        para[halfkeys[i]].append(key)
                        p.append(key)
                #    print(key)
        helperbool2=1
        i=i+1

    #----------------------------finding orthologs---------------------------

    ortho={}
    o=[]
    i=0
    cur="woop"
    helperbool=0
    helperbool2=0
    boollist=[]
    for x in range(amountofspec):
        boollist.append(0)
    for key in keys:

        if cur != halfkeys[i]:
            cur=halfkeys[i]
            #helperbool=0
            ortho[halfkeys[i]]=[]
            for x in range(amountofspec):
                ortho[halfkeys[i]].append("NA")
                boollist[x]=0
            #print(key)
            #print(tax2[key])
            #print(taxorder[tax2[key]])
            boollist[taxorder[tax2[key]]]=1


        if boollist[taxorder[tax2[key]]]==0:
            if float(fkeysE[key]) <= evalue:
                #print(key)
                if key not in o:
                    ortho[halfkeys[i]][taxorder[tax2[key]]]=key
                    o.append(key)

            #    print(tax2[key])
            #    print(taxorder[tax2[key]])
                boollist[taxorder[tax2[key]]]=1

        i=i+1





    #--------------------------finding co-ortholongs ---------------------------

    coortho={}
    for key in indhalfkeys:
        coortho[key]=[]
        if key in ortho.keys():
            for okey in ortho[key]:
                if okey!="NA":
                    if seckeys[okey] in para.keys():
                        for pkey in para[seckeys[okey]]:
                            coortho[key].append(pkey)
                            #print(pkey)
    #-------------------spitting out the output---------------------------------











    helperstring="fullrelationships.txt"

    output=open(helperstring, "w+")
    for key in indhalfkeys:
        output.write("Seq ID and File:   ")
        output.write(key)
        output.write("\n")
        output.write("PARALOGS----------------------------------------------")
        output.write("\n")
        for pkey in para[key]:
            output.write(seckeys[pkey])
            output.write("\t")
            output.write(fkeysE[(str(key)+"\t"+str(seckeys[pkey]))])
            output.write("\n")
        output.write("ORTHOLOGS---------------------------------------------")
        output.write("\n")
        for okey in ortho[key]:
            if okey!="NA":
                output.write(seckeys[okey])
                output.write("\t")
                output.write(fkeysE[(str(key)+"\t"+str(seckeys[okey]))])
                output.write("\n")
        output.write("CO-ORTHOLOGS------------------------------------------")
        output.write("\n")
        for ckey in coortho[key]:
            if ckey!="NA":
                output.write(seckeys[ckey])
                output.write("\t")
                t=(str(key)+"\t"+str(seckeys[ckey]))
                if t in fkeysE.keys():
                    output.write(fkeysE[t])
                output.write("\n")
        output.write("\n")
        output.write("\n")

    output.close()



















    loglist=[0,0]
    indexlist=[0,0]
    for key in indhalfkeys:
        for orth in ortho[key]:
            if orth!="NA":
                #print(orth)
                #print(keysE[orth])
                if float(fkeysE[orth]) != 0.0:
                    e=-(math.log10(float(fkeysE[orth])))

                    loglist[0]=loglist[0]+e
                else:
                    loglist[0]=loglist[0]+200
                    e=200
                indexlist[0]=indexlist[0]+1
                #print(e)
        for corth in coortho[key]:

            #print(corth)
            if float(fkeysE[corth]) != 0.0:
                e=-(math.log10(float(fkeysE[corth])))
            #print(e)
                loglist[0]=loglist[0]+e
            else:
                loglist[0]=loglist[0]+200
                e=200
            indexlist[0]=indexlist[0]+1
        for par in para[key]:

            #print(corth)
            if float(fkeysE[par]) != 0.0:
                e=-(math.log10(float(fkeysE[par])))
            #print(e)
                loglist[1]=loglist[1]+e
            else:
                loglist[1]=loglist[1]+200
                e=200
            indexlist[1]=indexlist[1]+1
    #print(loglist[0])
    #print(indexlist[0])
    normortho=loglist[0]/indexlist[0]
    normpara=loglist[1]/indexlist[1]



    output=open("graph.abc", "w+")
    fw={}

    network={}
    for key in indhalfkeys:
        vert=[]
        for orth in ortho[key]:
            if orth!="NA":
                if float(fkeysE[orth]) != 0.0:
                    e=-(math.log10(float(fkeysE[orth])))
                    e=e/normortho
                else:
                    e=200/normortho

                woop=str(orth+'\t'+str(e))
                output.write(woop)
                output.write("\n")
                edge=Edgy(orth,e)
                vert.append(edge)
                fw[orth]=e
                #print(edge.getL())
                #print(edge.getW())
        for corth in coortho[key]:
            if float(fkeysE[corth]) != 0.0:
                e=-(math.log10(float(fkeysE[corth])))
                e=e/normortho
            else:
                e=200/normortho
            woop=str(corth+'\t'+str(e))
            output.write(woop)
            output.write("\n")
            edge=Edgy(corth,e)
            vert.append(edge)
            fw[corth]=e
            #print(edge.getL())
            #print(edge.getW())
        for par in para[key]:
            if float(fkeysE[par]) != 0.0:
                e=-(math.log10(float(fkeysE[par])))
                e=e/normpara
            else:
                e=200/normpara
            woop=str(par+'\t'+str(e))
            output.write(woop)
            output.write("\n")
            edge=Edgy(par,e)
            vert.append(edge)
            fw[par]=e
            #print(edge.getL())
            #print(edge.getW())

        network[key]=vert

    output.close()

    os.system("mcxload -abc graph.abc -o OUTPUT.mci -write-tab OUTPUT.tab")

    com = "mcl OUTPUT.mci -I "
    com = com + str(inflat)
    com = com + " -o OUTPUT.out"

    os.system(com)
    #--------------------

    #---------------------
    tab2prot={}
    keykeys=[]

    helperfile=open("OUTPUT.tab","r")
    lines = helperfile.readlines()
    for line in lines:
        prot=""
        tab=""
        helperbool=0
        for char in line:
            if char == '\t':
                helperbool=1
                continue
            if helperbool==0:
                tab=tab+char
            else:
                if char=='\n':
                    continue
                prot=prot+char
        #print(tab)
        #print(prot)
        tab2prot[tab]=prot
        keykeys.append(prot)


    prot2cluster={}

    id2prot={}
    prota=[]
    helperfile=open("OUTPUT.tab","r")
    lines = helperfile.readlines()
    for line in lines:
        prot=""
        id=""
        helperbool=0
        b2=0
        b3=0
        for char in line:
            if char == '\t':
                helperbool=1
                continue


            if helperbool==1:

                if char=='\n':
                    continue
                if b2==1 and char =='|':
                    b3=1
                if b2==1 and b3==0:
                    id=id+char
                if char=='|':
                    b2=1
                prot=prot+char
        #print(tab)
        #print(prot)
        #print(id)

        id2prot[id]=prot
        prota.append(prot)


    helperfile.close()

    helperfile=open("OUTPUT.out","r")
    lines = helperfile.readlines()
    x=0
    prot=""
    cluster=""
    prot2clus={}
    for line in lines:
        if x<7:
            x=x+1
            continue
        prot=""
        helperbool=0
        hb=0
        if list(line)[0]!=' ':
            #prot=""
            cluster=""
            #helperbool=0
            #hb=0
            if list(line)[0] ==')':
                break
            for char in line:
                if char==" " or char=='\n':
                    if helperbool==1:
                        if hb==1:
                            if prot != "$":
                                if tab2prot[prot] not in prot2clus.keys():
                                    print(prot)
                                    prot2clus[tab2prot[prot]]=[]
                                prot2clus[tab2prot[prot]].append(cluster)
                                #print(tab2prot[prot])
                                #print(cluster)
                    prot=""
                    helperbool=1
                    continue
                if helperbool==0:
                    cluster=cluster+char
                else:
                    hb=1
                    prot=prot+char
        else:
            x=0
            prot=""
            for char in line:

                if x<7:
                    x=x+1
                    continue

                if char==" " or char=='\n':
                    if prot != "$":
                        if tab2prot[prot] not in prot2clus.keys():
                            print(prot)
                            prot2clus[tab2prot[prot]]=[]
                        prot2clus[tab2prot[prot]].append(cluster)
                        #print(tab2prot[prot])
                        #print(cluster)
                        prot=""
                else:
                    prot=prot+char

    goodpara=[]
    goodorth=[]
    goodco=[]
    i1=0
    i2=0
    for key in indhalfkeys:

        for pkey in para[key]:
            if firkeys[pkey]!=seckeys[pkey]:
                i1=i1+1
                if firkeys[pkey] in prot2clus.keys() and seckeys[pkey] in prot2clus.keys():
                    for x in prot2clus[firkeys[pkey]]:
                        for y in prot2clus[seckeys[pkey]]:
                            if x == y:
    #                        print(x)
    #                        print(y)
    #                        print(pkey)
                                goodpara.append(pkey)
                                i2=i2+1
        for ckey in coortho[key]:
            if firkeys[ckey]!=seckeys[ckey]:
                i1=i1+1
                if firkeys[ckey] in prot2clus.keys() and seckeys[ckey] in prot2clus.keys():
                    for x in prot2clus[firkeys[ckey]]:
                        for y in prot2clus[seckeys[ckey]]:
                            if x == y:
                                goodco.append(ckey)
                                i2=i2+1


        for okey in ortho[key]:
            if okey != "NA" :
                i1=i1+1
                if firkeys[okey] in prot2clus.keys() and seckeys[okey] in prot2clus.keys():
                    for x in prot2clus[firkeys[okey]]:
                        for y in prot2clus[seckeys[okey]]:
                            if x == y:
                                goodorth.append(okey)
                                i2=i2+1



    print(i1)
    print(i2)




    helperfile=open("txtHelper.txt","r")
    lines = helperfile.readlines()
    fullrec=[]
    strchecker="woop"
    send={}
    sstart={}
    for line in lines:
        b1=0
        b2=0
        b3=0
        b4=0
        se=""
        ss=""
        fid=""
        string=""
        for char in line:
            if char !='\t' and b1==0:
                fid=fid+char
            if b1==1:
                if char !='\t' and b2==0:
                    string=string+char
                if b2 == 1:

                    if char != '\t'and b3==0:
                        ss=ss+char
                    if b3==1:
                        if char != '\t':
                            se=se+char
                        else:
                            if string not in sstart.keys():
                                sstart[string]=[]
                                send[string]=[]
                            sstart[string].append(ss)
                            send[string].append(se)
                            if strchecker!=fid:
                                strchecker=fid
                                fullrec.append(string)

                            #print(string)
                            #print(ss)
                            #print(se)
                            break
                    if char == '\t':
                        b3=1
                if char == '\t':
                    b2=1
            if char == '\t':
                b1 =1

    #print(fullrec)




    goodclus=[]



    ##this one ------------------------------------------------------------------

    for key in fquary:
        if key in id2prot.keys():
            if id2prot[key] in prot2clus.keys():
                for clus in prot2clus[id2prot[key]]:
                    if clus not in goodclus:

                        goodclus.append(clus)

    #print(goodclus)
    forth=[]
    fpara=[]
    fco=[]
    f1=0
    f2=0
    f3=0
    for orth in goodorth:
        for cluster in goodclus:
            for y in prot2clus[seckeys[orth]]:
                if y==cluster:
                    forth.append(orth)
                    f1=f1+1
    for para in goodpara:
        for cluster in goodclus:
            for y in prot2clus[seckeys[para]]:
                if y==cluster:
                    fpara.append(para)
                    f2=f2+1
    for co in goodco:
        for cluster in goodclus:
            for y in prot2clus[seckeys[co]]:
                if y==cluster:
                    fco.append(co)
                    f3=f3+1
    print(len(forth)+len(fpara)+len(fco))


    idfir={}
    idsec={}
    megalistf=[]
    megalistf.extend(forth)
    megalistf.extend(fpara)
    megalistf.extend(fco)




    for key in megalistf:
        pre=""
        fir=""
        sec=""
        ib=""
        b1=0
        b2=0
        b3=0
        for char in key:
            if char !='|' and b1==0:
                pre=pre+char
            if b1==1:
                if char !='|' and b2==0:
                    fir=fir+char
                if b2 == 1:

                    if char != '|'and b3==0:
                        ib=ib+char
                    if b3==1:
                        if char != '|':
                            sec=sec+char
                        else:
                            idfir[key]=fir
                            idsec[key]=sec

                            #print(string)
                            #print(ss)
                            #print(se)
                            break
                    if char == '|':
                        b3=1
                if char == '|':
                    b2=1
            if char == '|':
                b1 =1


    #for key in indhalfkeys:
    #    print(key)
    #    print(sstart[key])
    #    print()












    helperstring= "preresultF.txt"
    output=open(helperstring, "w+")
    output.write("id_1\tid_2\tE-valueF\tWeight-F\tRelationship-F\t")
    output.write("\n")
    for key in megalistf:

        output.write(idfir[key])
        output.write("\t")




        output.write(idsec[key])
        output.write("\t")


        output.write(str(fkeysE[key]))
        output.write("\t")


        output.write(str(fw[key]))
        output.write("\t")
        b=0
        if key in fco:
            output.write("Co-ortholog")
            b=1
        if key in fpara and b==0:
            output.write("In-paralog")
            b=1
        if key in forth and b==0:
            output.write("Ortholog")
        output.write("\n")



    output.close()





















































































if domseq !="":



    #-------------------------domain code---------------------------------------

    homology_search(domseq,"firstDB","domtxtHelper.txt")


    helperfile=open("domtxtHelper.txt","r")
    lines = helperfile.readlines()
    fullrec=[]
    strchecker="woop"
    send={}
    sstart={}
    domain={}
    for line in lines:
        b1=0
        b2=0
        b3=0
        b4=0
        se=""
        ss=""
        fid=""
        string=""
        for char in line:
            if char !='\t' and b1==0:
                fid=fid+char
            if b1==1:
                if char !='\t' and b2==0:
                    string=string+char
                if b2 == 1:

                    if char != '\t'and b3==0:
                        ss=ss+char
                    if b3==1:
                        if char != '\t':
                            se=se+char
                        else:
                            if string not in sstart.keys():
                                sstart[string]=[]
                                send[string]=[]
                            sstart[string].append(int(ss))
                            send[string].append(int(se))
                            if strchecker!=fid:
                                strchecker=fid
                                fullrec.append(string)
                                domain[string]=fid
                                #print(fid)
                                #print(string)
                            #print(string)
                            #print(ss)
                            #print(se)
                            break
                    if char == '\t':
                        b3=1
                if char == '\t':
                    b2=1
            if char == '\t':
                b1 =1
    #print(fullrec)

    sseqids=[]
    tax={}
    fulltax={}



    curfull="woop"

    for line in lines:
        b=0
        b2=0
        fcur=""
        fkey=""
        ftax=""
        for char in line:
            if char=='\t':
                fkey=fcur
            if b == 0:
                if char != '\t':
                    fcur=fcur+char
            else:
                if char=='#':
                    b2=1

                if b2==1:
                    if char=='\t':
                        break
                    if char != '\t':
                        ftax=ftax+char
            if char=='\t':
                b=1
        if curfull!=fkey:
            curfull=fkey
            fulltax[fkey]=ftax





    index=0

    for line in lines:  #                for each line in the result in the blastp
        l=list(line)
        index=0
        tempstring=""
        helperbool=0
        b=0
        for chara in l:  #                    for each chararicter in the line
            if chara == '\t':
                for i in range(len(l)):#        record teh line until the tab char
                    if l[index+i+1] != '\t':
                        tempstring= tempstring + l[index+i+1]
                    else:
                        for sseqid in sseqids: #                   for each sseqid
                            if str(sseqid) == tempstring:
                                helperbool=1
                        if helperbool ==0:
                            sseqids.append(tempstring)#add ids to the array of ids
                            #print(tempstring)
                        break
                break
            index=index+1


    dba={}
    seq={}
    fullrec=[]

    for record in rec:
        string=(record.seq)
        id=record.id
        dba[id]=[]
        seq[id]=string
        for char in string:
            dba[id].append(False)


    for key in sstart.keys():
        index=0
        for y in sstart[key]:
            for x in range((int(sstart[key][index])-1),((send[key][index])-1)):
                dba[key][x]=True
            index=index+1
    newrec=[]
    strchecker="woop"

    for key in dba.keys():


        dc=0
        dom=""
        x=0
        bool=0
        index1=1
        index2=0
        for char in dba[key]:
            index1=index1+1
            if char==True:
                index2=index2+1
                if bool==0:
                    dc=dc+1
                dom=dom+seq[key][x]
                if len(dba[key])<=x+1:
                    continue
                if dba[key][x+1]==False:
                    name="dom"+str(dc)+'-'+str(index1-index2)+'-'+str(index1)+"*"+str(key)
                    #print(name)
                    fullrec.append(name)
                    record=SeqRecord(Seq(dom),id=name)
                    newrec.append(record)
                    index2=0
                bool=1
            else:
                bool=0
                dom=""
            x=x+1







    '''

    for key in sstart.keys():
        start=(2**31-1.)
        i=i+1
        for x in sstart[key]:
            if int(x) < start:
                start=int(x)
        sstart[key]=start
        end=0
        for x in send[key]:
            if int(x) > end:
                end=int(x)
        send[key]=end
        #print(key)
        #print(sstart[key])
        #print(send[key])
    print(i)


    newrec=[]
    for record in rec:


        for sseqid in sseqids:
            if record.id==sseqid:
                #print("-------------------------")
                #print(record.seq)
                record=SeqRecord(Seq(str(record.seq)[(sstart[sseqid]-1):(send[sseqid]+1)]), id=sseqid)
                #print(record.id)
                #print(str(record.seq)[(sstart[sseqid]-1):(send[sseqid]+1)])
                #print("---")
                #print(record.seq)
                #print("-------------------------")
                newrec.append(record)
     '''

    print(len(rec))
    print(len(sseqids))
    print(len(newrec))



    #print(newrec)

    #----------------read and writes the file for the second blastP---------------
    newref=open("homologyD.fasta","w+")
    SeqIO.write(newrec,"homologyD.fasta","fasta")
    newref.close()

    build_blast_db("homologyD.fasta","homologyD")

    #----------------blastp of all vs all----------------------------------------
    print("allvsall start")


    homology_search("homologyD.fasta","homologyD","all_vs_allD.blastp")


    print("allvsall finish")





    helperfile=open("all_vs_allD.blastp","r")
    lines = helperfile.readlines()

    keys=[]
    halfkeys=[]
    indhalfkeys=[]
    dkeysE={}
    keysB={}
    keysP={}
    minE={}
    maxB={}
    maxP={}
    tax={}
    seckeys={}
    firkeys={}
    samebool=[]
    taxorder={}
    #sstart={}
    #send={}
    taxorderhelper=0

    #aveP={}
    #-------------go through the files and find all the info of the -------------
    #--------------------AllvsAll Blastp outputfile------------------------------

    #          things you should know about about blastp files
    #          there differnt info is sesperated by a tab (\t)
    #       in the code when you see tab(\t) its the end thing inbetween info


    for line in lines:# for each line
        l=list(line)
        #print(l)
        #helper varibles
        index=0
        fin =0
        f=""
        s=""
        firststring=""
        tempstring=""
        helperbool=0
        taxbool=0
        taxa=""

        for chara in l: #                      for each character in each line
            if chara !='\t':#     if not the end of the line add to temp sting
                tempstring= tempstring +chara
                if helperbool == 0:#          find first string and add to var
                    firststring= firststring+ chara
                    f=f+chara
                else:
                    s=s+chara
                #print(tempstring)
            if chara == '\t':#   if end of string in file add taxonomy to dic
                taxbool=0
                tax[tempstring]=taxa
                if not taxa in taxorder.keys():#       checks for no reapeats
                    taxorder[taxa]=taxorderhelper
                    taxorderhelper=taxorderhelper+1
                if helperbool==1:#               if the second id was passed
                    if f == s: #      if first id and sencond id is the same
                        samebool.append(1)#to be honest i dont even use this
                    else: #         but to scared to scared to get rid of it
                        samebool.append(0)
                    index=index+1       #index helps counts the where we are
                    keys.append(tempstring) #    add the full key to the dic
                    charbool=0              #       We are done with id part
                    secstring=""
                    for char in tempstring: #   for the charartors in stirng
                        if charbool==1: #    if on second id place it in dic
                            secstring=secstring+char
                        if char=='\t':
                            charbool=1
                    firkeys[tempstring]=firststring
                    seckeys[tempstring]=secstring

                    if not firststring in indhalfkeys:#checking for repeats
                        indhalfkeys.append(firststring)
    #                      make a array of half keys to help later funtions
    #                                       Both unquie and non unquie list
                    halfkeys.append(firststring)
                    #if not tempstring in aveB.keys():
                    #    aveB[tempstring]=[]
                    #    aveP[tempstring]=[]

            #------ after the ID and into the other data------
                    #ss=""
                    #if secstring not in sstart.keys():
                    #    sstart[secstring]=[]
                    for i in range(len(l)):# the start info and redoring it
                        if l[index]!='\t':
                    #        ss=ss+l[index]
                            index=index+1
                        if l[index]=='\t':
                            index=index+1
                            break
                    #if ss not in sstart[secstring]:
                    #sstart[secstring].append(ss)
                    #se=""
                    #if secstring not in send.keys():
                    #    send[secstring]=[]
                    for i in range(len(l)):#  the end info and recording it
                        if l[index]!='\t':
                    #        se=se+l[index]
                            index=index+1
                        if l[index]=='\t':
                            index=index+1
                            break
                    #if se not in send[secstring]:
                    #send[secstring].append(se)
                    e=""#                      the E-value and recording it
                    for i in range(len(l)):
                        if l[index]!='\t':
                            e=e+l[index]
                            index=index+1
                        if l[index]=='\t':
                            dkeysE[tempstring]=e
                            if not firststring in minE.keys():
                                minE[firststring]=e
                            else:
                                if float(e) < float(minE[firststring]):
                                    minE[firststring]=e
                            index=index+1
                            break
                    b=""#                   the bit-score and recording it
                    for i in range(len(l)):
                        if l[index]!='\t':
                            b=b+l[index]
                            index=index+1
                        if l[index]=='\t':
                            keysB[tempstring]=b
                            if not firststring in maxB.keys():
                                maxB[firststring]=b
                            else:
                                if float(b) > float(maxB[firststring]):
                                    maxB[firststring]=b

                            index=index+1
                            break
                    p=""#                    The postive and recording it
                    for i in range(len(l)):#   for the length of the line
                        if l[index]!='\t':#            tab inbetween info
                            p=p+l[index]
                            index=index+1 #             add more to index
                        if index==len(l):#                 if end of line
                            keysP[tempstring]=p#
                            if not firststring in maxP.keys():
                                maxP[firststring]=p#       add pos to dic
                            else:
                                if float(p) > float(maxP[firststring]):
                                    maxP[firststring]=p
                            fin=1
                            # helper bool that tells us that we are done
                            index=index+1 #          gotta add the index
                            break
                tempstring=tempstring+'\t'
                helperbool=1
                taxa=taxa+'#'

            if fin == 1:#                                  if done break
                break

            if taxbool==1:     #                   stuff for the tax dic
                taxa=taxa+chara
            if chara=='#':
                taxbool =1
            index=index+1

    helperfile.close()


    #print(len(keys))
    #for key in keys:
    #    print(key)

    print("allvall scanner all done")




    tax2={}#                               helps dics for the taxonomy
    tax1={}
    taxorder={}
    taxorderhelper=0

    for key in keys: #                                for the all keys
        t1=""#                                               first tax
        t2=""#                                              second tax
        taxa=tax[key]#      gets the tax dic and slits it into two dic
        helperbool=0
        for char in taxa:
            if char=='#':
                helperbool=1
                continue
            if helperbool == 0:
                t1=t1+char
            if helperbool ==1:
                t2=t2+char
        if not t2 in taxorder.keys():
            taxorder[t2]=taxorderhelper
            taxorderhelper=taxorderhelper+1
        tax1[key]=t1
        tax2[key]=t2
    dtax1=tax1
    dtax2=tax2
    #-----------------------finding paralogs------------------------------------
    para={}
    i=0
    cur="woop"
    helperbool=0
    helperbool2=0
    for key in keys:

        if cur != halfkeys[i]:
            cur=halfkeys[i]
            helperbool=0
            para[halfkeys[i]]=[]
            helperbool2=0


        if tax2[key]!=tax1[key]:
            helperbool=1
        if helperbool == 0:
            if helperbool2==1:
                #print(key)
                #print(keysE[key])
                if float(dkeysE[key]) <= evalue:
                    para[halfkeys[i]].append(key)
                #    print(key)
        helperbool2=1
        i=i+1

    #----------------------------finding orthologs---------------------------

    ortho={}
    i=0
    cur="woop"
    helperbool=0
    helperbool2=0
    boollist=[]
    for x in range(amountofspec):
        boollist.append(0)
    for key in keys:

        if cur != halfkeys[i]:
            cur=halfkeys[i]
            #helperbool=0
            ortho[halfkeys[i]]=[]
            for x in range(amountofspec):
                ortho[halfkeys[i]].append("NA")
                boollist[x]=0
            #print(key)
            #print(tax2[key])
            #print(taxorder[tax2[key]])
            boollist[taxorder[tax2[key]]]=1


        if boollist[taxorder[tax2[key]]]==0:
            if float(dkeysE[key]) <= evalue:
                #print(key)
                ortho[halfkeys[i]][taxorder[tax2[key]]]=key

            #    print(tax2[key])
            #    print(taxorder[tax2[key]])
                boollist[taxorder[tax2[key]]]=1

        i=i+1


    #--------------------------finding co-ortholongs ---------------------------

    coortho={}

    for key in indhalfkeys:
        coortho[key]=[]
        for okey in ortho[key]:
            if okey!="NA":
                if seckeys[okey] in para.keys():
                    for pkey in para[seckeys[okey]]:

                        coortho[key].append(pkey)
                        #print(pkey)
    #-------------------spitting out the output---------------------------------










    helperstring="domainrelationships.txt"
    output=open(helperstring, "w+")
    for key in indhalfkeys:
        output.write("Seq ID and File:   ")
        output.write(key)
        output.write("\n")
        output.write("PARALOGS----------------------------------------------")
        output.write("\n")
        for pkey in para[key]:
            output.write(seckeys[pkey])
            output.write("\t")
            output.write(dkeysE[(str(key)+"\t"+str(seckeys[pkey]))])
            output.write("\n")
        output.write("ORTHOLOGS---------------------------------------------")
        output.write("\n")
        for okey in ortho[key]:
            if okey!="NA":
                output.write(seckeys[okey])
                output.write("\t")
                output.write(dkeysE[(str(key)+"\t"+str(seckeys[okey]))])
                output.write("\n")
        output.write("CO-ORTHOLOGS------------------------------------------")
        output.write("\n")
        for ckey in coortho[key]:
            if ckey!="NA":
                output.write(seckeys[ckey])
                output.write("\t")
                t=(str(key)+"\t"+str(seckeys[ckey]))
                if t in dkeysE.keys():
                    output.write(dkeysE[t])
                output.write("\n")
        output.write("\n")
        output.write("\n")

    output.close()








    loglist=[0,0]
    indexlist=[0,0]
    for key in indhalfkeys:
        for orth in ortho[key]:
            if orth!="NA":
                #print(orth)
                #print(keysE[orth])
                if float(dkeysE[orth]) != 0.0:
                    e=-(math.log10(float(dkeysE[orth])))

                    loglist[0]=loglist[0]+e
                else:
                    loglist[0]=loglist[0]+200
                    e=200
                indexlist[0]=indexlist[0]+1
                #print(e)
        for corth in coortho[key]:

            #print(corth)
            if float(dkeysE[corth]) != 0.0:
                e=-(math.log10(float(dkeysE[corth])))
            #print(e)
                loglist[0]=loglist[0]+e
            else:
                loglist[0]=loglist[0]+200
                e=200
            indexlist[0]=indexlist[0]+1
        for par in para[key]:

            #print(corth)
            if float(dkeysE[par]) != 0.0:
                e=-(math.log10(float(dkeysE[par])))
            #print(e)
                loglist[1]=loglist[1]+e
            else:
                loglist[1]=loglist[1]+200
                e=200
            indexlist[1]=indexlist[1]+1
    #print(loglist[0])
    #print(indexlist[0])
    normortho=loglist[0]/indexlist[0]
    normpara=loglist[1]/indexlist[1]

    output=open("graphD.abc", "w+")

    dw={}
    network={}
    for key in indhalfkeys:
        vert=[]
        for orth in ortho[key]:
            if orth!="NA":
                if float(dkeysE[orth]) != 0.0:
                    e=-(math.log10(float(dkeysE[orth])))
                    e=e/normortho
                else:
                    e=200/normortho

                woop=str(orth+'\t'+str(e))
                output.write(woop)
                output.write("\n")
                edge=Edgy(orth,e)
                vert.append(edge)
                dw[orth]=e
                #print(edge.getL())
                #print(edge.getW())
        for corth in coortho[key]:
            if float(dkeysE[corth]) != 0.0:
                e=-(math.log10(float(dkeysE[corth])))
                e=e/normortho
            else:
                e=200/normortho
            woop=str(corth+'\t'+str(e))
            output.write(woop)
            output.write("\n")
            edge=Edgy(corth,e)
            vert.append(edge)
            dw[corth]=e
            #print(edge.getL())
            #print(edge.getW())
        for par in para[key]:
            if float(dkeysE[par]) != 0.0:
                e=-(math.log10(float(dkeysE[par])))
                e=e/normpara
            else:
                e=200/normpara
            woop=str(par+'\t'+str(e))
            output.write(woop)
            output.write("\n")
            edge=Edgy(par,e)
            vert.append(edge)
            dw[par]=e
            #print(edge.getL())
            #print(edge.getW())

        network[key]=vert

    output.close()

    os.system("mcxload -abc graphD.abc -o OUTPUTd.mci -write-tab OUTPUTd.tab")

    com = "mcl OUTPUTd.mci -I "
    com = com + str(inflat)
    com = com + " -o OUTPUTd.out"

    os.system(com)


    #---------------------
    tab2prot={}
    keykeys=[]

    helperfile=open("OUTPUTd.tab","r")
    lines = helperfile.readlines()
    for line in lines:
        prot=""
        tab=""
        helperbool=0
        for char in line:
            if char == '\t':
                helperbool=1
                continue
            if helperbool==0:
                tab=tab+char
            else:
                if char=='\n':
                    continue
                prot=prot+char
        #print(tab)
        #print(prot)
        #print(tab)
        tab2prot[tab]=prot
        keykeys.append(prot)



    id2prot={}
    prota=[]
    helperfile=open("OUTPUTd.tab","r")
    lines = helperfile.readlines()
    for line in lines:
        prot=""
        id=""
        helperbool=0
        b2=0
        b3=0
        for char in line:
            if char == '\t':
                helperbool=1
                continue


            if helperbool==1:

                if char=='\n':
                    continue
                if b2==1 and char =='|':
                    b3=1
                if b2==1 and b3==0:
                    id=id+char
                if char=='|':
                    b2=1
                prot=prot+char
        #print(tab)
        #print(prot)
        #print(id)

        if id not in id2prot.keys():
            id2prot[id]=[]
        id2prot[id].append(prot)
        prota.append(prot)





    helperfile.close()

    helperfile=open("OUTPUTd.out","r")
    lines = helperfile.readlines()
    x=0
    prot=""
    cluster=""
    prot2clus={}
    for line in lines:
        if x<7:
            x=x+1
            continue
        prot=""
        helperbool=0
        hb=0
        if list(line)[0]!=' ':
            #prot=""
            cluster=""
            #helperbool=0
            #hb=0
            if list(line)[0] ==')':
                break
            for char in line:
                if char==" " or char=='\n':
                    if helperbool==1:
                        if hb==1:
                            if prot != "$":
                                if tab2prot[prot] not in prot2clus.keys():
                                    print(prot)
                                    prot2clus[tab2prot[prot]]=[]
                                prot2clus[tab2prot[prot]].append(cluster)
                                #print(tab2prot[prot])
                                #print(cluster)
                    prot=""
                    helperbool=1
                    continue
                if helperbool==0:
                    cluster=cluster+char
                else:
                    hb=1
                    prot=prot+char
        else:
            x=0
            prot=""
            for char in line:

                if x<7:
                    x=x+1
                    continue

                if char==" " or char=='\n':
                    if prot != "$":
                        if tab2prot[prot] not in prot2clus.keys():
                            print(prot)
                            prot2clus[tab2prot[prot]]=[]
                            #print(tab2prot[prot])
                        prot2clus[tab2prot[prot]].append(cluster)
                        #print(tab2prot[prot])
                        #print(cluster)
                        prot=""
                else:
                    prot=prot+char

    goodpara=[]
    goodorth=[]
    goodco=[]
    i1=0
    i2=0
    for key in indhalfkeys:

        for pkey in para[key]:
            if firkeys[pkey]!=seckeys[pkey]:
                i1=i1+1
                if firkeys[pkey] in prot2clus.keys() and seckeys[pkey] in prot2clus.keys():
                    for x in prot2clus[firkeys[pkey]]:
                        for y in prot2clus[seckeys[pkey]]:
                            if x == y:
    #                        print(x)
    #                        print(y)
    #                        print(pkey)
                                goodpara.append(pkey)
                                i2=i2+1
        for ckey in coortho[key]:
            if firkeys[ckey]!=seckeys[ckey]:
                i1=i1+1
                if firkeys[ckey] in prot2clus.keys() and seckeys[ckey] in prot2clus.keys():
                    for x in prot2clus[firkeys[ckey]]:
                        for y in prot2clus[seckeys[ckey]]:
                            if x == y:
                                goodco.append(ckey)
                                i2=i2+1


        for okey in ortho[key]:
            if okey != "NA" :
                i1=i1+1
                if firkeys[okey] in prot2clus.keys() and seckeys[okey] in prot2clus.keys():
                    for x in prot2clus[firkeys[okey]]:
                        for y in prot2clus[seckeys[okey]]:
                            if x == y:
                                goodorth.append(okey)
                                i2=i2+1


    print(len(goodpara))
    print(len(goodorth))
    print(len(goodco))

    goodclus=[]

    for key in dquary:
        if key in id2prot.keys():
            for prot in id2prot[key]:
                if prot in prot2clus.keys():
                    for clus in prot2clus[prot]:
                        if clus not in goodclus:

                            goodclus.append(clus)

    print(len(goodclus))

    #print(goodclus)
    dorth=[]
    dpara=[]
    dco=[]
    d1=0
    d2=0
    d3=0
    for orth in goodorth:
        for cluster in goodclus:
            for y in prot2clus[seckeys[orth]]:
                if y==cluster:
                    dorth.append(orth)
                    d1=d1+1
    for para in goodpara:
        for cluster in goodclus:
            for y in prot2clus[seckeys[para]]:
                if y==cluster:
                    dpara.append(para)
                    d2=d2+1
    for co in goodco:
        for cluster in goodclus:
            for y in prot2clus[seckeys[co]]:
                if y==cluster:
                    dco.append(co)
                    d3=d3+1


    megalistd=[]
    megalistd.extend(dorth)
    megalistd.extend(dpara)
    megalistd.extend(dco)
    idfir={}
    idsec={}
    print(len(megalistd))

    for key in megalistd:
        pre=""
        fir=""
        sec=""
        ib=""
        b1=0
        b2=0
        b3=0
        for char in key:
            if char !='|' and b1==0:
                pre=pre+char
            if b1==1:
                if char !='|' and b2==0:
                    fir=fir+char
                if b2 == 1:

                    if char != '|'and b3==0:
                        ib=ib+char
                    if b3==1:
                        if char != '|':
                            sec=sec+char
                        else:
                            idfir[key]=fir
                            idsec[key]=sec

                            #print(string)
                            #print(ss)
                            #print(se)
                            break
                    if char == '|':
                        b3=1
                if char == '|':
                    b2=1
            if char == '|':
                b1 =1


    helperstring= "preresultD.txt"
    output=open(helperstring, "w+")
    output.write("id_1\tdomain1\tid_2\tdomain2\tE-valueD\tWeight-D\tRelationship-D\t")
    output.write("\n")
    for key in megalistd:

        output.write(idfir[key])
        output.write("\t")
        if firkeys[key] in domain.keys():
            output.write(domain[firkeys[key]])
        else:
            hs=""
            hs=hs+idfir[key]
            hs=hs+"_domain1"

            output.write(hs)
        output.write("\t")
        output.write(idsec[key])
        output.write("\t")
        if seckeys[key] in domain.keys():
            output.write(domain[seckeys[key]])
        else:
            hs=""
            hs=hs+idsec[key]
            hs=hs+"_domain1"

            output.write(hs)
        output.write("\t")

        output.write(str(dkeysE[key]))
        output.write("\t")
        output.write(str(dw[key]))
        output.write("\t")

        b=0
        if key in dco:
            output.write("Co-ortholog")
            b=1
        if key in dpara and b==0:
            output.write("In-paralog")
            b=1
        if key in dorth and b==0:
            output.write("Ortholog")

        output.write("\n")

    output.close()






ds=0

if domseq != "":
    ds=1
    prot1={}
    prot2={}

    finallist=[]
    dom1={}
    dom2={}
    doms={}
    bases1={}
    bases2={}
    i=0
    for key in megalistd:

        modkey=""
        fdom=""
        sdom=""
        base1=""
        mod1=""
        base2=""
        mod2=""
        b1=0
        b2=0
        b3=0

        for char in key:
            if char !='*' and b1==0 :
                mod1=mod1+char
                fdom=fdom+char
            if b1==1:
                if char !='\t' and b2==0:
                    base1=base1+char
                    mod1=mod1+char
                    modkey=modkey+char
                if b2 == 1:

                    if char != '*'and b3==0:
                        mod2=mod2+char
                        sdom=sdom+char

                    if b3==1:
                        modkey=modkey+char
                        mod2=mod2+char
                        base2=base2+char

                    if char == '*':
                        if b3==0:
                            mod2=mod2+'*'
                        b3=1

                if char == '\t':
                    b2=1
                    modkey=modkey+'\t'
                    #base=base+'\t'
                    #mod=mod+'\t'
            if char == '*':
                if b1==0:
                    mod1=mod1+'*'
                b1 =1


        #print(modkey)
        #print(mod2)
        #print(mod1)
        doms[key]=modkey
        dom1[key]=fdom
        dom2[key]=sdom
        bases1[mod1]=base1
        bases2[mod2]=base2
        i=i+1
        p1=""
        p2=""
        b4=0
        for char in base1:
            if char == '#':
                b4=1

            if b4==0:
                p1=p1+char
        b4=0
        for char in base2:
            if char == '#':
                b4=1
            if b4==0:
                p2=p2+char
        prot1[key]=p1
        prot2[key]=p2
        if fullseq != "":
            if modkey in megalistf:

                finallist.append(key)


        else:
            finallist.append(key)






fs=0
if domseq == "":
    prot1={}
    prot2={}
    finallist=[]

if fullseq != "":
    fs=1

    for key in megalistf:
        b1=0
        b2=0
        b3=0
        p1=""
        p2=""
        for char in key:
            if char !='#' and b1==0:
                p1=p1+char
            if b1==1:

                if b2 == 1:

                    if char != '#'and b3==0:
                        p2=p2+char
                    if b3==1:
                        b3=1
                    if char == '#':
                        b3=1
                if char == '\t':
                    b2=1
            if char == '#':
                b1 =1

        prot1[key]=p1
        prot2[key]=p2

        if domseq=="":
            finallist.append(key)



#print(finallist)
#print(i)

woop = set()
for key in finallist:
    key2=seckeys[key]+'\t'+firkeys[key]
    if key not in woop and key2 not in woop:
        woop.add(key)

finallist=list(woop)
if fullseq != "":
    woop = set()
    for key in megalistf:
        key2=seckeysf[key]+'\t'+firkeysf[key]
        if key not in woop and key2 not in woop:
            woop.add(key)

    megalistf=list(woop)



if domseq != "":
    woop = set()
    for key in megalistd:
        key2=seckeys[key]+'\t'+firkeys[key]
        if key not in woop and key2 not in woop:
            woop.add(key)

    megalistd=list(woop)


if out=="":
    output=open("results.txt", "w+")
else:
    helperstring= out + ".txt"
    output=open(helperstring, "w+")



#-------------------------------------------------------------------------
#non domain and full orthologs






















































































if domseq == ""  and fullseq =="":
    print("all vs all start (this might take a hot sec)")
    homology_search("ref.fasta","firstDB","all_vs_allnone.blastp")



    print("all vs all done ")



    print("reading the all vs all output")





    helperfile=open("all_vs_allnone.blastp","r")
    lines = helperfile.readlines()

    keys=[]
    halfkeys=[]
    indhalfkeys=[]
    fkeysE={}
    keysB={}
    keysP={}
    minE={}
    maxB={}
    maxP={}
    tax={}
    seckeys={}
    firkeys={}
    samebool=[]
    taxorder={}
    last="woop"

    #sstart={}
    #send={}
    taxorderhelper=0

    #aveP={}
    #-------------go through the files and find all the info of the -------------
    #--------------------AllvsAll Blastp outputfile------------------------------

    #          things you should know about about blastp files
    #          there differnt info is sesperated by a tab (\t)
    #       in the code when you see tab(\t) its the end thing inbetween info


    for line in lines:# for each line
        l=list(line)
        #print(l)
        #helper varibles
        index=0
        fin =0
        f=""
        s=""
        firststring=""
        tempstring=""
        helperbool=0
        taxbool=0
        taxa=""

        for chara in l: #                      for each character in each line
            if chara !='\t':#     if not the end of the line add to temp sting
                tempstring= tempstring +chara
                if helperbool == 0:#          find first string and add to var
                    firststring= firststring+ chara
                    f=f+chara
                else:
                    s=s+chara
                #print(tempstring)
            if chara == '\t':#   if end of string in file add taxonomy to dic
                taxbool=0
                tax[tempstring]=taxa
                if not taxa in taxorder.keys():#       checks for no reapeats
                    taxorder[taxa]=taxorderhelper
                    taxorderhelper=taxorderhelper+1
                if helperbool==1:#               if the second id was passed
                    if f == s: #      if first id and sencond id is the same
                        samebool.append(1)#to be honest i dont even use this
                    else: #         but to scared to scared to get rid of it
                        samebool.append(0)
                    index=index+1       #index helps counts the where we are
                    keys.append(tempstring) #    add the full key to the dic

                    charbool=0              #       We are done with id part
                    secstring=""
                    for char in tempstring: #   for the charartors in stirng
                        if charbool==1: #    if on second id place it in dic
                            secstring=secstring+char
                        if char=='\t':
                            charbool=1
                    firkeys[tempstring]=firststring
                    seckeys[tempstring]=secstring

                    if not firststring in indhalfkeys:#checking for repeats
                        indhalfkeys.append(firststring)
    #                      make a array of half keys to help later funtions
    #                                       Both unquie and non unquie list
                    halfkeys.append(firststring)
                    #if not tempstring in aveB.keys():
                    #    aveB[tempstring]=[]
                    #    aveP[tempstring]=[]

            #------ after the ID and into the other data------
                    #ss=""
                    #if secstring not in sstart.keys():
                    #    sstart[secstring]=[]
                    for i in range(len(l)):# the start info and redoring it
                        if l[index]!='\t':
                    #        ss=ss+l[index]
                            index=index+1
                        if l[index]=='\t':
                            index=index+1
                            break
                    #if ss not in sstart[secstring]:
                    #sstart[secstring].append(ss)
                    #se=""
                    #if secstring not in send.keys():
                    #    send[secstring]=[]
                    for i in range(len(l)):#  the end info and recording it
                        if l[index]!='\t':
                    #        se=se+l[index]
                            index=index+1
                        if l[index]=='\t':
                            index=index+1
                            break
                    #if se not in send[secstring]:
                    #send[secstring].append(se)
                    e=""#                      the E-value and recording it
                    for i in range(len(l)):
                        if l[index]!='\t':
                            e=e+l[index]
                            index=index+1
                        if l[index]=='\t':
                            fkeysE[tempstring]=e
                            if not firststring in minE.keys():
                                minE[firststring]=e
                            else:
                                if float(e) < float(minE[firststring]):
                                    minE[firststring]=e
                            index=index+1
                            break
                    b=""#                   the bit-score and recording it
                    for i in range(len(l)):
                        if l[index]!='\t':
                            b=b+l[index]
                            index=index+1
                        if l[index]=='\t':
                            keysB[tempstring]=b
                            if not firststring in maxB.keys():
                                maxB[firststring]=b
                            else:
                                if float(b) > float(maxB[firststring]):
                                    maxB[firststring]=b

                            index=index+1
                            break
                    p=""#                    The postive and recording it
                    for i in range(len(l)):#   for the length of the line
                        if l[index]!='\t':#            tab inbetween info
                            p=p+l[index]
                            index=index+1 #             add more to index
                        if index==len(l):#                 if end of line
                            keysP[tempstring]=p#
                            if not firststring in maxP.keys():
                                maxP[firststring]=p#       add pos to dic
                            else:
                                if float(p) > float(maxP[firststring]):
                                    maxP[firststring]=p
                            fin=1
                            # helper bool that tells us that we are done
                            index=index+1 #          gotta add the index
                            break
                tempstring=tempstring+'\t'
                helperbool=1
                taxa=taxa+'#'

            if fin == 1:#                                  if done break
                break

            if taxbool==1:     #                   stuff for the tax dic
                taxa=taxa+chara
            if chara=='#':
                taxbool =1
            index=index+1

    helperfile.close()




    print("allvall scanner all done")
    #----------------------creats helpers comp sturctures-----------------------

    tax2={}#                               helps dics for the taxonomy
    tax1={}
    taxorder={}
    taxorderhelper=0

    for key in keys: #                                for the all keys
        t1=""#                                               first tax
        t2=""#                                              second tax
        taxa=tax[key]#      gets the tax dic and slits it into two dic
        helperbool=0
        for char in taxa:
            if char=='#':
                helperbool=1
                continue
            if helperbool == 0:
                t1=t1+char
            if helperbool ==1:
                t2=t2+char
        if not t2 in taxorder.keys():
            taxorder[t2]=taxorderhelper
            taxorderhelper=taxorderhelper+1
        tax1[key]=t1
        tax2[key]=t2

    ftax1=tax1
    ftax2=tax2

    #-----------------------finding paralogs------------------------------------
    para={}
    p=[]
    i=0
    cur="woop"
    helperbool=0
    helperbool2=0
    for key in keys:

        if cur != halfkeys[i]:
            cur=halfkeys[i]
            helperbool=0
            para[halfkeys[i]]=[]
            helperbool2=0


        if tax2[key]!=tax1[key]:
            helperbool=1
        if helperbool == 0:
            if helperbool2==1:
                #print(key)
                #print(keysE[key])
                if float(fkeysE[key]) <= evalue:
                    if key not in p:
                        para[halfkeys[i]].append(key)
                        p.append(key)
                #    print(key)
        helperbool2=1
        i=i+1

    #----------------------------finding orthologs---------------------------

    ortho={}
    o=[]
    i=0
    cur="woop"
    helperbool=0
    helperbool2=0
    boollist=[]
    for x in range(amountofspec):
        boollist.append(0)
    for key in keys:

        if cur != halfkeys[i]:
            cur=halfkeys[i]
            #helperbool=0
            ortho[halfkeys[i]]=[]
            for x in range(amountofspec):
                ortho[halfkeys[i]].append("NA")
                boollist[x]=0
            #print(key)
            #print(tax2[key])
            #print(taxorder[tax2[key]])
            boollist[taxorder[tax2[key]]]=1


        if boollist[taxorder[tax2[key]]]==0:
            if float(fkeysE[key]) <= evalue:
                #print(key)
                if key not in o:
                    ortho[halfkeys[i]][taxorder[tax2[key]]]=key
                    o.append(key)

            #    print(tax2[key])
            #    print(taxorder[tax2[key]])
                boollist[taxorder[tax2[key]]]=1

        i=i+1





    #--------------------------finding co-ortholongs ---------------------------

    coortho={}
    for key in indhalfkeys:
        coortho[key]=[]
        if key in ortho.keys():
            for okey in ortho[key]:
                if okey!="NA":
                    if seckeys[okey] in para.keys():
                        for pkey in para[seckeys[okey]]:
                            coortho[key].append(pkey)
                            #print(pkey)
    #-------------------spitting out the output---------------------------------











    helperstring="fullrelationships.txt"

    output=open(helperstring, "w+")
    for key in indhalfkeys:
        output.write("Seq ID and File:   ")
        output.write(key)
        output.write("\n")
        output.write("PARALOGS----------------------------------------------")
        output.write("\n")
        for pkey in para[key]:
            output.write(seckeys[pkey])
            output.write("\t")
            output.write(fkeysE[(str(key)+"\t"+str(seckeys[pkey]))])
            output.write("\n")
        output.write("ORTHOLOGS---------------------------------------------")
        output.write("\n")
        for okey in ortho[key]:
            if okey!="NA":
                output.write(seckeys[okey])
                output.write("\t")
                output.write(fkeysE[(str(key)+"\t"+str(seckeys[okey]))])
                output.write("\n")
        output.write("CO-ORTHOLOGS------------------------------------------")
        output.write("\n")
        for ckey in coortho[key]:
            if ckey!="NA":
                output.write(seckeys[ckey])
                output.write("\t")
                t=(str(key)+"\t"+str(seckeys[ckey]))
                if t in fkeysE.keys():
                    output.write(fkeysE[t])
                output.write("\n")
        output.write("\n")
        output.write("\n")

    output.close()



















    loglist=[0,0]
    indexlist=[0,0]
    for key in indhalfkeys:
        for orth in ortho[key]:
            if orth!="NA":
                #print(orth)
                #print(keysE[orth])
                if float(fkeysE[orth]) != 0.0:
                    e=-(math.log10(float(fkeysE[orth])))

                    loglist[0]=loglist[0]+e
                else:
                    loglist[0]=loglist[0]+200
                    e=200
                indexlist[0]=indexlist[0]+1
                #print(e)
        for corth in coortho[key]:

            #print(corth)
            if float(fkeysE[corth]) != 0.0:
                e=-(math.log10(float(fkeysE[corth])))
            #print(e)
                loglist[0]=loglist[0]+e
            else:
                loglist[0]=loglist[0]+200
                e=200
            indexlist[0]=indexlist[0]+1
        for par in para[key]:

            #print(corth)
            if float(fkeysE[par]) != 0.0:
                e=-(math.log10(float(fkeysE[par])))
            #print(e)
                loglist[1]=loglist[1]+e
            else:
                loglist[1]=loglist[1]+200
                e=200
            indexlist[1]=indexlist[1]+1
    #print(loglist[0])
    #print(indexlist[0])
    normortho=loglist[0]/indexlist[0]
    normpara=loglist[1]/indexlist[1]



    output=open("graphn.abc", "w+")
    fw={}

    network={}
    for key in indhalfkeys:
        vert=[]
        for orth in ortho[key]:
            if orth!="NA":
                if float(fkeysE[orth]) != 0.0:
                    e=-(math.log10(float(fkeysE[orth])))
                    e=e/normortho
                else:
                    e=200/normortho

                woop=str(orth+'\t'+str(e))
                output.write(woop)
                output.write("\n")
                edge=Edgy(orth,e)
                vert.append(edge)
                fw[orth]=e
                #print(edge.getL())
                #print(edge.getW())
        for corth in coortho[key]:
            if float(fkeysE[corth]) != 0.0:
                e=-(math.log10(float(fkeysE[corth])))
                e=e/normortho
            else:
                e=200/normortho
            woop=str(corth+'\t'+str(e))
            output.write(woop)
            output.write("\n")
            edge=Edgy(corth,e)
            vert.append(edge)
            fw[corth]=e
            #print(edge.getL())
            #print(edge.getW())
        for par in para[key]:
            if float(fkeysE[par]) != 0.0:
                e=-(math.log10(float(fkeysE[par])))
                e=e/normpara
            else:
                e=200/normpara
            woop=str(par+'\t'+str(e))
            output.write(woop)
            output.write("\n")
            edge=Edgy(par,e)
            vert.append(edge)
            fw[par]=e
            #print(edge.getL())
            #print(edge.getW())

        network[key]=vert

    output.close()

    os.system("mcxload -abc graphn.abc -o OUTPUTN.mci -write-tab OUTPUTN.tab")

    com = "mcl OUTPUTN.mci -I "
    com = com + str(inflat)
    com = com + " -o OUTPUTN.out"

    os.system(com)
    #--------------------

    #---------------------
    tab2prot={}
    keykeys=[]

    helperfile=open("OUTPUTN.tab","r")
    lines = helperfile.readlines()
    for line in lines:
        prot=""
        tab=""
        helperbool=0
        for char in line:
            if char == '\t':
                helperbool=1
                continue
            if helperbool==0:
                tab=tab+char
            else:
                if char=='\n':
                    continue
                prot=prot+char
        #print(tab)
        #print(prot)
        tab2prot[tab]=prot
        keykeys.append(prot)


    prot2cluster={}

    id2prot={}
    prota=[]
    helperfile=open("OUTPUTN.tab","r")
    lines = helperfile.readlines()
    for line in lines:
        prot=""
        id=""
        helperbool=0
        b2=0
        b3=0
        for char in line:
            if char == '\t':
                helperbool=1
                continue


            if helperbool==1:

                if char=='\n':
                    continue
                if b2==1 and char =='|':
                    b3=1
                if b2==1 and b3==0:
                    id=id+char
                if char=='|':
                    b2=1
                prot=prot+char
        #print(tab)
        #print(prot)
        #print(id)

        id2prot[id]=prot
        prota.append(prot)


    helperfile.close()

    helperfile=open("OUTPUTN.out","r")
    lines = helperfile.readlines()
    x=0
    prot=""
    cluster=""
    prot2clus={}
    for line in lines:
        if x<6:
            x=x+1
            continue
        prot=""
        helperbool=0
        hb=0
        if list(line)[0]!=' ':
            #prot=""
            cluster=""
            #helperbool=0
            #hb=0
            if list(line)[0] ==')':
                break
            for char in line:
                if char==" " or char=='\n':
                    if helperbool==1:
                        if hb==1:
                            if prot != "$":
                                if tab2prot[prot] not in prot2clus.keys():
                                    print(prot)
                                    prot2clus[tab2prot[prot]]=[]
                                prot2clus[tab2prot[prot]].append(cluster)
                                #print(tab2prot[prot])
                                #print(cluster)
                    prot=""
                    helperbool=1
                    continue
                if helperbool==0:
                    cluster=cluster+char
                else:
                    hb=1
                    prot=prot+char
        else:
            x=0
            prot=""
            for char in line:

                if x<6:
                    x=x+1
                    continue

                if char==" " or char=='\n':
                    if prot != "$":

                        if tab2prot[prot] not in prot2clus.keys():
                            print(prot)
                            prot2clus[tab2prot[prot]]=[]
                        prot2clus[tab2prot[prot]].append(cluster)
                        #print(tab2prot[prot])
                        #print(cluster)
                        prot=""
                else:
                    prot=prot+char

    goodpara=[]
    goodorth=[]
    goodco=[]
    i1=0
    i2=0
    for key in indhalfkeys:

        for pkey in para[key]:
            if firkeys[pkey]!=seckeys[pkey]:
                i1=i1+1
                if firkeys[pkey] in prot2clus.keys() and seckeys[pkey] in prot2clus.keys():
                    for x in prot2clus[firkeys[pkey]]:
                        for y in prot2clus[seckeys[pkey]]:
                            if x == y:
    #                        print(x)
    #                        print(y)
    #                        print(pkey)
                                goodpara.append(pkey)
                                i2=i2+1
        for ckey in coortho[key]:
            if firkeys[ckey]!=seckeys[ckey]:
                i1=i1+1
                if firkeys[ckey] in prot2clus.keys() and seckeys[ckey] in prot2clus.keys():
                    for x in prot2clus[firkeys[ckey]]:
                        for y in prot2clus[seckeys[ckey]]:
                            if x == y:
                                goodco.append(ckey)
                                i2=i2+1


        for okey in ortho[key]:
            if okey != "NA" :
                i1=i1+1
                if firkeys[okey] in prot2clus.keys() and seckeys[okey] in prot2clus.keys():
                    for x in prot2clus[firkeys[okey]]:
                        for y in prot2clus[seckeys[okey]]:
                            if x == y:
                                goodorth.append(okey)
                                i2=i2+1







    '''
    helperfile=open("txtHelperN.txt","r")
    lines = helperfile.readlines()
    fullrec=[]
    strchecker="woop"
    send={}
    sstart={}
    for line in lines:
        b1=0
        b2=0
        b3=0
        b4=0
        se=""
        ss=""
        fid=""
        string=""
        for char in line:
            if char !='\t' and b1==0:
                fid=fid+char
            if b1==1:
                if char !='\t' and b2==0:
                    string=string+char
                if b2 == 1:

                    if char != '\t'and b3==0:
                        ss=ss+char
                    if b3==1:
                        if char != '\t':
                            se=se+char
                        else:
                            if string not in sstart.keys():
                                sstart[string]=[]
                                send[string]=[]
                            sstart[string].append(ss)
                            send[string].append(se)
                            if strchecker!=fid:
                                strchecker=fid
                                fullrec.append(string)

                            #print(string)
                            #print(ss)
                            #print(se)
                            break
                    if char == '\t':
                        b3=1
                if char == '\t':
                    b2=1
            if char == '\t':
                b1 =1

    #print(fullrec)





    forth=[]
    fpara=[]
    fco=[]
    f1=0
    f2=0
    f3=0
    for orth in goodorth:
        for cluster in goodclus:
            for y in prot2clus[seckeys[orth]]:
                if y==cluster:
                    forth.append(orth)
                    f1=f1+1
    for para in goodpara:
        for cluster in goodclus:
            for y in prot2clus[seckeys[para]]:
                if y==cluster:
                    fpara.append(para)
                    f2=f2+1
    for co in goodco:
        for cluster in goodclus:
            for y in prot2clus[seckeys[co]]:
                if y==cluster:
                    fco.append(co)
                    f3=f3+1
    print(len(forth)+len(fpara)+len(fco))'''














    fco=goodco
    fpara=goodpara
    forth=goodorth










    idfir={}
    idsec={}
    megalistn=[]
    megalistn.extend(goodorth)
    megalistn.extend(goodpara)
    megalistn.extend(goodco)




    for key in megalistn:
        pre=""
        fir=""
        sec=""
        ib=""
        b1=0
        b2=0
        b3=0
        for char in key:
            if char !='|' and b1==0:
                pre=pre+char
            if b1==1:
                if char !='|' and b2==0:
                    fir=fir+char
                if b2 == 1:

                    if char != '|'and b3==0:
                        ib=ib+char
                    if b3==1:
                        if char != '|':
                            sec=sec+char
                        else:
                            idfir[key]=fir
                            idsec[key]=sec

                            #print(string)
                            #print(ss)
                            #print(se)
                            break
                    if char == '|':
                        b3=1
                if char == '|':
                    b2=1
            if char == '|':
                b1 =1


    #for key in indhalfkeys:
    #    print(key)
    #    print(sstart[key])
    #    print()












    helperstring= "preresultF.txt"
    output=open(helperstring, "w+")
    output.write("id_1\tid_2\tE-valueF\tWeight-F\tRelationship-F\t")
    output.write("\n")
    for key in megalistn:

        output.write(idfir[key])
        output.write("\t")




        output.write(idsec[key])
        output.write("\t")


        output.write(str(fkeysE[key]))
        output.write("\t")


        output.write(str(fw[key]))
        output.write("\t")
        b=0
        if key in fco:
            output.write("Co-ortholog")
            b=1
        if key in fpara and b==0:
            output.write("In-paralog")
            b=1
        if key in forth and b==0:
            output.write("Ortholog")
        output.write("\n")



    output.close()



    for key in megalistn:
        b1=0
        b2=0
        b3=0
        p1=""
        p2=""
        for char in key:
            if char !='#' and b1==0:
                p1=p1+char
            if b1==1:

                if b2 == 1:

                    if char != '#'and b3==0:
                        p2=p2+char
                    if b3==1:
                        b3=1
                    if char == '#':
                        b3=1
                if char == '\t':
                    b2=1
            if char == '#':
                b1 =1

        prot1[key]=p1
        prot2[key]=p2







































if out=="":
    helperstring="results.txt"
else:
    helperstring=out+".txt"
output=open(helperstring, "w+")





if fullseq != "" and domseq !="":


    output.write("Species_1\tProtein_1\tDomain_1\tSpecies_2\tProtein_2\tDomain_2\tE-value_Full\t E-value_Domain\tWeight_Full\tWeight_Domain\tRelationship_Full\tRelationship_Domain")
    output.write("\n")
    output.write("------------------------------------------------------------------------------------------------------------------------------------------------------------------")
    output.write("\n")
    for key in finallist:
        output.write(str(ftax1[doms[key]]))
        output.write("\t")
        output.write(prot1[key])
        output.write("\t")
        hs=dom1[key]

        output.write(hs)
        output.write("\t")
        output.write(str(ftax2[doms[key]]))
        output.write("\t")
        output.write(prot2[key])
        output.write("\t")
        hs=dom2[key]

        output.write(hs)
        output.write("\t")
        output.write(str(fkeysE[doms[key]]))
        output.write("\t")
        output.write(str(dkeysE[key]))
        output.write("\t")
        output.write(str(fw[doms[key]]))
        output.write("\t")
        output.write(str(dw[key]))
        output.write("\t")
        b=0
        if doms[key] in fpara and ftax1[doms[key]]==ftax2[doms[key]]:
            output.write("In-paralog")
            b=1
        if doms[key] in fco and b==0 and ftax1[doms[key]]!=ftax2[doms[key]]:
            output.write("Co-ortholog")
            b=1
        if doms[key] in forth and b==0:
            output.write("Ortholog")
        output.write("\t")
        b=0
        if key in dpara and dtax1[key]==dtax2[key]:
            output.write("In-paralog")
            b=1
        if key in dco and b==0 and dtax1[key]!=dtax2[key]:
            output.write("Co-ortholog")
            b=1
        if key in dorth and b==0:
            output.write("Ortholog")

        output.write("\n")

if fullseq != "":


    output.write("------------------------------------------------------------------------------------------------------------------------------------------------------------------")
    output.write("\n")
    output.write("FULL GRAPH")
    output.write("\n")
    output.write("Species_1\tProtein_1\tSpecies_2\tProtein_2\tE-value_Full\tWeight_Full\tRelationship_Full")
    output.write("\n")
    output.write("------------------------------------------------------------------------------------------------------------------------------------------------------------------")
    output.write("\n")
    for key in megalistf:
        output.write(str(ftax1[key]))
        output.write("\t")
        output.write(prot1[key])
        output.write("\t")
        output.write(str(ftax2[key]))
        output.write("\t")
        output.write(prot2[key])
        output.write("\t")
        output.write(str(fkeysE[key]))
        output.write("\t")
        output.write(str(fw[key]))
        output.write("\t")
        b=0
        if key in fpara and ftax1[key]==ftax2[key]:
            output.write("In-paralog")
            b=1
        if key in fco and b==0 and ftax1[key]!=ftax2[key]:
            output.write("Co-ortholog")
            b=1
        if key in forth and b==0:
            output.write("Ortholog")
        output.write("\t")
        output.write("\n")

if domseq != "":
    output.write("------------------------------------------------------------------------------------------------------------------------------------------------------------------")
    output.write("\n")
    output.write("Domain GRAPH")
    output.write("\n")
    output.write("Species_1\tProtein_1\tDomain_1\tSpecies_2\tProtein_2\tDomain_2\tE-value_Domain\tWeight_Domain\tRelationship_Domain")
    output.write("\n")
    output.write("------------------------------------------------------------------------------------------------------------------------------------------------------------------")
    output.write("\n")
    for key in megalistd:
        output.write(str(dtax1[key]))
        output.write("\t")
        output.write(prot1[key])
        output.write("\t")

        hs=dom1[key]

        output.write(hs)
        output.write("\t")
        output.write(str(dtax2[key]))
        output.write("\t")
        output.write(prot2[key])
        output.write("\t")
        hs=dom2[key]

        output.write(hs)
        output.write("\t")
        output.write(str(dkeysE[key]))
        output.write("\t")
        output.write(str(dw[key]))
        output.write("\t")
        b=0
        if key in dpara and dtax1[key]==dtax2[key]:
            output.write("In-paralog")
            b=1
        if key in dco and b==0 and dtax1[key]!=dtax2[key]:
            output.write("Co-ortholog")
            b=1
        if key in dorth and b==0:
            output.write("Ortholog")

        output.write("\n")






if fullseq == "" and domseq == "" :


    output.write("------------------------------------------------------------------------------------------------------------------------------------------------------------------")
    output.write("\n")
    output.write("None quary  GRAPH")
    output.write("\n")
    output.write("Species_1\tProtein_1\tSpecies_2\tProtein_2\tE-value_Full\tWeight_Full\tRelationship_Full")
    output.write("\n")
    output.write("------------------------------------------------------------------------------------------------------------------------------------------------------------------")
    output.write("\n")
    for key in megalistn:
        output.write(str(ftax1[key]))
        output.write("\t")
        output.write(prot1[key])
        output.write("\t")
        output.write(str(ftax2[key]))
        output.write("\t")
        output.write(prot2[key])
        output.write("\t")
        output.write(str(fkeysE[key]))
        output.write("\t")
        output.write(str(fw[key]))
        output.write("\t")
        b=0
        if key in fpara and ftax1[key]==ftax2[key]:
            output.write("In-paralog")
            b=1
        if key in fco and b==0 and ftax1[key]!=ftax2[key]:
            output.write("Co-ortholog")
            b=1
        if key in forth and b==0:
            output.write("Ortholog")
        output.write("\t")
        output.write("\n")













output.close()










#for key in finallist:
#    print(key)
#    print(fkeysE[key])
#    print(dkeysE[key])
#    print(fw[key])
#    print(dw[key])
#    if firkeys[key] in domain.keys():
#        print(domain[firkeys[key]])
#    if seckeys[key] in domain.keys():
#        print(domain[seckeys[key]])




#for key in megalistf:
#    print(key)
#    print(idfir[key])
#    print(idsec[key])
#print("-----------------------------------------------------------------------")
#for key in megalistd:
#    print(key)
#    print(idfir[key])
#    print(idsec[key])
#    print(domain[firkeys[key]])
#    print(domain[seckeys[key]])







'''
print(dorth)
print(dpara)
print(dco)
print(forth)
print(fpara)
print(fco)
print(f1)
print(f2)
print(f3)
print(d1)
print(d2)
print(d3)
'''
'''

start=(2**31-1.)
i=0
for key in indhalfkeys:
    start=(2**31-1.)
    i=i+1
    for x in sstart[key]:
        if int(x) < start:
            start=int(x)
    sstart[key]=start
    end=0
    for x in send[key]:
        if int(x) > end:
            end=int(x)
    send[key]=end
    print(key)
    print(sstart[key])
    print(send[key])
print(i)


'''



'''
tempstr=""
templist=[]

for key in fullrec:
    tempstr=""
    tempstr=key+fulltax[key]
    templist.append(tempstr)
fullrec=templist
#print(fullrec)
'''
#records = list(SeqIO.parse(goodfile, "fasta"))
#i2=0
#for x in goodpara:
    #print(x)
#    i2=i2+1
#for x in goodorth:
    #print(x)
#    i2=i2+1
#for x in goodco:
    #print(x)
#    i2=i2+1

#$print(i1)
#print(i2)
'''
if out=="":
    output=open("results.txt", "w+")
else:
    helperstring= out + ".txt"
    output=open(helperstring, "w+")
for key in indhalfkeys:
    output.write("Seq ID and File:   ")
    output.write(key)
    output.write("\n")
    output.write("PARALOGS----------------------------------------------")
    output.write("\n")
    for pkey in para[key]:
        output.write(seckeys[pkey])
        output.write("\n")
    output.write("ORTHOLOGS---------------------------------------------")
    output.write("\n")
    for okey in ortho[key]:
        if okey!="NA":
            output.write(seckeys[okey])
            output.write("\n")
    output.write("CO-ORTHOLOGS------------------------------------------")
    output.write("\n")
    for ckey in coortho[key]:
        if ckey!="NA":
            output.write(seckeys[ckey])
            output.write("\n")
    output.write("\n")
    output.write("\n")

output.close()
print(para)
#-(math.log10())

'''

'''unqcombo=[]
indexdic={}
normalize={}
for key in keys:
    if tax[key] in unqcombo:
        indexdic[tax[key]]=indexdic[tax[key]]+1
        normalize[tax[key]]=normalize[tax[key]]-(math.log10(keysE[key]))
    else:
        unqcombo.append(tax[key])
        indexdic[tax[key]]=1
        normalize[tax[key]]=(-(math.log10(keysE[key])))
'''

import os
import pandas as pd

seqsero2output = snakemake.input[0]
forward = snakemake.input[1]
reverse = snakemake.input[2]
styphi = snakemake.config["readthreshold"]["styphi"]
biphasic = snakemake.config["readthreshold"]["biphasic"]
threads = snakemake.threads


#save the file names as forward and reverse

serotype = ["4:i:1,2","I 4,[5],12:i:-"]

#create name of the input/output file
splittedpath = forward.split("/")
name = splittedpath[-1].replace('_R1.fastq.gz','')

with open(seqsero2output) as results:
    serotypedata = results.read()
    if any(x in serotypedata for x in serotype):
        
        #create amplicons list to loop over
        amplicons = [snakemake.input[3], snakemake.input[4]]
        
        #for loop for both amplicons
        for amplicon in amplicons:
            ampliconnamesplit = amplicon.split("/")
            ampliconname = ampliconnamesplit[-1].replace('.fasta','')
            
            #run all the commands on commandline to generate depth file
            os.system("bwa index %s"%(amplicon))
            os.system("bwa mem -k 17 -t %s %s %s %s > %s.sam"%(threads,amplicon,forward,reverse,name))
            os.system("samtools sort -@ %s -n -O sam %s.sam | samtools fixmate -m -O bam - %s.bam"%(threads,name,name))
            os.system("rm %s.sam"%(name))
            os.system("samtools sort -@ %s -O bam -o %s.sorted.bam %s.bam"%(threads,name,name))
            os.system("rm %s.bam"%(name))
            os.system("samtools depth %s.sorted.bam | gzip > %s.%s.depth.txt.gz"%(name,name,ampliconname))
            os.system("rm %s.sorted.bam"%(name))
            os.system("gunzip %s.%s.depth.txt.gz"%(name,ampliconname))
        
            # open depth file as x
            with open('%s.%s.depth.txt' % (name, ampliconname)) as x:
                # if loop which depends on the amplicon
                if ampliconname == 'sense_59_antisense_83':
                    count = 0
                    # loop over every line in the depth file
                    for read in x:
                        # split the line based on tabs
                        splitread = read.split('\t')
                        # remove enter (\n)
                        splitread[2].replace('\n', '')
        
                        # if sense_59_antisense_83 amplicon is not present in WGS there is a gap between the positions 425 and 1130
                        if int(splitread[1]) > 425 and int(splitread[1]) < 1130:
                            count += int(splitread[2])
                            countsense_59_antisense_83 = count
                    # cutoff threshold for amount of reads between these positions is set to default 5000
                    if count > biphasic:
                        sense_59_antisense_83output = True
                    else:
                        sense_59_antisense_83output = False
                    # remove the depth file
                    os.system("rm %s.%s.depth.txt" % (name, ampliconname))
        
                # repeat the same as before but if it is for the FFLIB_FFLIA amplicon
                elif ampliconname == 'FFLIB_FFLIA':
                    count = 0
                    for read in x:
                        splitread = read.split('\t')
                        splitread[2].replace('\n', '')
                        # if FFLIB_FFLIA is not present in WGS there is a gap between the positions 425 and 1130
                        if int(splitread[1]) > 75 and int(splitread[1]) < 780:
                            count += int(splitread[2])
                            countFFLIB_FFLIA = count
                                    
                    if count > styphi:
                        FFLIB_FFLIAoutput = True
                    else:
                        FFLIB_FFLIAoutput = False
                    os.system("rm %s.%s.depth.txt" % (name, ampliconname))
                
    else:
        FFLIB_FFLIAoutput = False
        sense_59_antisense_83output = False
        countFFLIB_FFLIA = "N/A"
        countsense_59_antisense_83 = "N/A"

#determine if sample is biphasic, monophasic, or none
if FFLIB_FFLIAoutput == True and sense_59_antisense_83output == True:
    variant = "Biphasic"
elif FFLIB_FFLIAoutput == True and sense_59_antisense_83output == False:
    variant = "Monophasic"
else:
    variant = "N/A"
    
#add the output into a dataframe using pandas

seqsero_report = pd.read_csv(seqsero2output, sep='\t')
seqsero_report = seqsero_report.iloc[:,[0,3,4,5,6,7,8,9,10]]
seqsero_report.insert(7,"Monophasic/Biphasic", variant)
seqsero_report.insert(10,"FFLIB_FFLIAcount", countFFLIB_FFLIA)
seqsero_report.insert(11,"sense_59_antisense_83count", countsense_59_antisense_83)
seqsero_report.to_csv(name+"_combinedresult.tsv",sep='\t', index=False, header=True)

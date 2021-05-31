# CheckAmpliconsInSTyphi
This repository contains the additional parts that need to be added to the Juno pipeline in order to discriminate S. Typhimurium from other H:i serovars and monophasic S. Typhimurium from biphasic S. typhimurium

### Additions to main pipeline
- Added script (checkamplicons.py) to discriminate S. Typhimurium from other H:i serovars and monophasic S. Typhimurium from biphasic S. Typhimurium based on an in silico based approach
  - Script first determines if specific serovar is present based on the output of seqsero2, if not present, script returns a none value for that sample.
  - The script maps the amplicon sequences of the primers FFLIB and RFLIA, and sense_59 and antisense_83 against the WGS data. With the mapped reads, samtools is used to generate a depth file. 
  - The depth file shows the amount of reads for each nucleotide of the amplicon sequence. The script then counts the amount of reads between nucleotide 75 and 780 for the FFLIB and RFLIA amplicon and nucleotide 425 and 1130 for the sense_59 and antisense_83 amplicon
  - A parameter "ReadThreshold" is then used to determine if the sample is monophasic or biphasic S. typhimurium.
  - The output of this is a combined tsv file of seqsero2 and if the sample is monophasic, biphasic, or none
  - For our project, we have added 2 additional columns which can be removed later from the script. This contains the counts of the amount of reads in the regions which are described above.
- Added rule (checkamplicons) that runs the script that was described above
  - input is the result of seqsero2, WGS data, and the amplicon sequences
  - additional parameters are available, such as the amount of threads and the readthreshold
  - can run with the same version of tools that are used for the seqsero2 rule

### Changes in main pipeline
- Changed multireport rule and script (seqsero2_multireport.py). 
    - Multireport rule now takes a combined tsv file of seqsero2 and the checkamplicons rule.
    - Script now takes all columns and rows in the .iloc function. This is because the combined tsv file already filters out the unneccasary columns and adds the checkamplicon result.

### Possible errors
- ran into a problem when running snakemake using more cores and at the same time more threads for seqsero2. This is probably a local problem on our server.

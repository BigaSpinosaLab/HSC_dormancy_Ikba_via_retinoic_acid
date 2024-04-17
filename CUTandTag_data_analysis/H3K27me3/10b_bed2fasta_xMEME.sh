
# NOTE: This execution requires a local installation of MEME Suite

# Define variables
# Location with BEDs stored: Execute script '10a_BED_Exclusive_Regions_MEME.R' 
BED="Ikba_HSC_ROSHANA/BED_Exclusive/" 

# Location to store corresponding fasta files
OUT=$BED

# Location of ref genome sequence: Download ref genome Ensembl Release 102 mm10 
REF="/Volumes/projectscomput/cancer/db_files/Genomes/Ensembl/mouse/mm10/release-102/Mus_musculus.GRCm38.dna.primary_assembly.fa" 


for FILE in $BED/*.bed
do
    NAME=${FILE%.bed}
    SAMPLE=$(basename $NAME)
    
    /opt/local/libexec/meme-5.5.2/bed2fasta $FILE $REF  > $OUT/$SAMPLE.fa
    
done
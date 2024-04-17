
# Define variables
# Location with FASTA (derived from BED) stored

FASTA="Ikba_HSC_ROSHANA/BED_Exclusive"

# Location to store MEME Suite results
OUT="Ikba_HSC_ROSHANA/MEMESuite_results" 

# Location of databases: This can be retreived from MEME Suite
DB1="/Users/mmaqueda/Documents/Software/MEMESuite/motif_databases/MOUSE/HOCOMOCOv11_core_MOUSE_mono_meme_format.meme"

for FILE in $FASTA/*.fa
do
    NAME=${FILE%.fa}
    SAMPLE=$(basename $NAME)

    meme-chip -o $OUT/$SAMPLE  -minw 5 -maxw 16 -db $DB1 -meme-nmotifs 10 -streme-totallength 4000000 -fimo-skip -spamo-skip    $FILE

done



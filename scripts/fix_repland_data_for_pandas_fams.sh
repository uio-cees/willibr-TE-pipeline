#	Purpose: trim parseRM_GetLandscape.pl output for import to pandas, get class info.


FILEPATH=$1

cd $FILEPATH

#----------------------------------START OF LOOP-----------------------------------------#

for i in $(ls $FILEPATH |grep "Perc.all" );
do

#	Need total number of lines to "tail out" the correct amount of lines
LENGTH=$(cat $i | wc -l)
TAIL=$(($LENGTH-3))

grep -v "(" $i > $i.nosimple

#	This is the body of the data
column -t $i.nosimple | tail -n $TAIL > $i.repland


#	Getting header, removing "[" and "]"
head -n 2  $i | tail -n1 | \
sed 's/\[//g' | sed 's/\]//g' | column -t \
	>  header

#	Getting length of repland_complete
REPLAND_LENGTH=$(cat $i.repland | wc -l)

#	Attempting to create column file for later merge
yes ${i%%.*} | head -n $REPLAND_LENGTH > ${i%%.*}.left

#	Pasting left alias column with the correct data
paste ${i%%.*}.left $i.repland > ${i%%.*}.RL.tab


#	Deleting files
rm  $i.repland ${i%%.*}.left $i.nosimple

done
#---------------------------------END OF LOOP--------------------------------------------#

echo "ALIAS" $(cat header) > header_new

cat *RL.tab > teleost.RepeatLandscapes.body

awk '{$NF=""; print $0}' teleost.RepeatLandscapes.body \
	> teleost.body

cat header_new teleost.body > teleost.RepeatLandscapes.csv

rm teleost.RepeatLandscapes.tab
rm *RL.tab
rm header
rm header_new
rm teleost.body
rm teleost.RepeatLandscapes.body

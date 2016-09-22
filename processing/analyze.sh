for i in polerr2016-1 polerr2016-2
do
   java -Xmx60G -jar mageri-1.0.0.jar -I $i.json --import-preset preset.xml -O $i-out/ 2>&1 | tee $i-log.txt
done
YEAR2016=hists/WWW2016_v1.2.2/test2016_v6NmO
YEAR2017=hists/Loose2017_v5.1.4/test2017_v6NmO
YEAR2018=hists/Loose2018_v5.1.4/test2018_v6NmO

YEARCOMB=hists/combineyears/test_v6NmO

if [ ! -d "$YEARCOMB" ]; then
    echo "making dir $YEARCOMB"
    mkdir -p "$YEARCOMB"
fi

find "$YEAR2016" -maxdepth 1 -type f -name '*.root' -printf '%f\n' | while read line; do
    if [ ! -f "$YEAR2016/$line" ]; then
        echo "File $YEAR2016/$line does not exist. continue."
        continue
    fi
    if [ ! -f "$YEAR2017/$line" ]; then
        echo "File $YEAR2017/$line does not exist. continue."
        continue
    fi
    if [ ! -f "$YEAR2017/$line" ]; then
        echo "File $YEAR2017/$line does not exist. continue."
        continue
    fi
    echo "hadd -f $YEARCOMB/$line $YEAR2016/$line $YEAR2017/$line $YEAR2017/$line"
    hadd -f $YEARCOMB/$line $YEAR2016/$line $YEAR2017/$line $YEAR2017/$line
done


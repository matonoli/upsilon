OUT="effiLM""$1"".txt"
echo "$OUT"
touch "$OUT"
rm "$OUT"
touch "$OUT"
FILE="upsE""$1""_0420.root"
echo "$FILE"
root -l -q nTupEff.C\(\""$FILE"\","$2",0,2,9,\"tmp.txt\"\)
cat tmp.txt >> "$OUT"
root -l -q nTupEff.C\(\""$FILE"\","$2",0,2,5,\"tmp.txt\"\)
cat tmp.txt >> "$OUT"
root -l -q nTupEff.C\(\""$FILE"\","$2",0,5,7,\"tmp.txt\"\)
cat tmp.txt >> "$OUT"
root -l -q nTupEff.C\(\""$FILE"\","$2",0,7,9,\"tmp.txt\"\)
cat tmp.txt >> "$OUT"
OUT="effiH""$1"".txt"
echo "$OUT"
touch "$OUT"
rm "$OUT"
touch "$OUT"
root -l -q nTupEff.C\(\""$FILE"\","$2",1,2,9,\"tmp.txt\"\)
cat tmp.txt >> "$OUT"
root -l -q nTupEff.C\(\""$FILE"\","$2",1,2,5,\"tmp.txt\"\)
cat tmp.txt >> "$OUT"
root -l -q nTupEff.C\(\""$FILE"\","$2",1,5,7,\"tmp.txt\"\)
cat tmp.txt >> "$OUT"
root -l -q nTupEff.C\(\""$FILE"\","$2",1,7,9,\"tmp.txt\"\)
cat tmp.txt >> "$OUT"

echo "HI"
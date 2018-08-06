OUT="yields.txt"
echo "$OUT"
touch "$OUT"
rm "$OUT"
touch "$OUT"
root4star -l -q FitUpsilonMassAll.C\(2,9\)
cat tmp.txt >> "$OUT"
root4star -l -q FitUpsilonMassAll.C\(2,5\)
cat tmp.txt >> "$OUT"
root4star -l -q FitUpsilonMassAll.C\(5,7\)
cat tmp.txt >> "$OUT"
root4star -l -q FitUpsilonMassAll.C\(7,9\)
cat tmp.txt >> "$OUT"

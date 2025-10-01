from collections import Counter

S="ACGGGCATATGCGC"
alphabet=set(S)
counts=Counter(S)
percentages={}
for k in counts:
    percentages[k]=round(counts[k]/len(S)*100,2)

print("Alphabet:", alphabet)
print("Percentages:", percentages)

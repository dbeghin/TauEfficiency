import csv

tripletList=[]
with open("different.csv") as ffile:
    reader=csv.reader(ffile, delimiter=',')
    for row in reader:
        tripletList.append( (int(row[0]), int(row[1]), int(row[2])) )

newTriplet = sorted(tripletList, key=lambda triplet: triplet[0])

file_out = open("different_sorted.csv", "w")
for i in range(0, len(tripletList)):
    file_out.write(str(newTriplet[i][0])+", "+str(newTriplet[i][1])+", "+str(newTriplet[i][2])+"\n")

file_out.close()

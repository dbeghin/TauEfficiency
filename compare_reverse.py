import csv

Andrew=[]
with open("Andrew.csv") as ffile:
    reader=csv.reader(ffile, delimiter=',')
    for row in reader:
        Andrew.append(row)


different=[]
file_out = open("different_reverse.csv", "w")
with open("Diego.csv") as ffile:
    reader=csv.reader(ffile, delimiter=',')
    for row in reader:
        if not row in Andrew:
            print row
            different.append(row)
            file_out.write(row[0]+", ")
            file_out.write(row[1]+", ")
            file_out.write(row[2])
            file_out.write("\n")

file_out.close()

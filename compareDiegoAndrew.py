import csv

Diego=[]
with open("Diego.csv") as ffile:
    reader=csv.reader(ffile, delimiter=',')
    for row in reader:
        Diego.append(row)


different=[]
file_out = open("different.csv", "w")
with open("Andrew.csv") as ffile:
    reader=csv.reader(ffile, delimiter=',')
    for row in reader:
        if not row in Diego:
            print row
            different.append(row)
            file_out.write(row[0]+", ")
            file_out.write(row[1]+", ")
            file_out.write(row[2])
            file_out.write("\n")

file_out.close()

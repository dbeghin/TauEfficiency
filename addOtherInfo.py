import csv

AndrewTotal=[]
with open("AndrewTotal.csv") as ffile:
    reader=csv.reader(ffile, delimiter=',')
    for row in reader:
        temp_list=[]
        for j in range(0, 9):
            temp_list.append(row[j])
        AndrewTotal.append(temp_list)


file_out = open("different_totalinfo.csv", "w")
with open("different_sorted.csv") as ffile:
    reader=csv.reader(ffile, delimiter=',')
    for row in reader:
        for i in range(0, len(AndrewTotal)):
            match = True
            for j in range(0, len(row)):
                if int(row[j]) != int(AndrewTotal[i][j]):
                    match = False
                    break
            if match:
                file_out.write(AndrewTotal[i][0])
                for j in range(1, len(AndrewTotal[i])):
                    file_out.write(","+AndrewTotal[i][j])
                file_out.write("\n")
                print AndrewTotal[i]
                break

file_out.close()

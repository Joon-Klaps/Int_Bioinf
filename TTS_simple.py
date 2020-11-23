path = "F://Desktop//bedFiles//vir_10.bed"

lists = []
listStartEnd = []
listEndTemp = []
file = open(path, "r")
lines = file.readlines()
listEndFinal = []

for line in lines:
    lists.append(line.split("\t"))

for i in range(0, len(lists)-1):
    listTemp = []
    listTemp.append(lists[i][1])
    listTemp.append(lists[i][2])
    listStartEnd.append(listTemp)
    listEndTemp.append(int(lists[i][2]))

print(listEndTemp[0])
print(type(listEndTemp[0]))

number = input("clustering close sites: ")
counts = input("Number of counts allowed: ")

listEndTemp.sort()

grouping = [[listEndTemp[0]]]
for i in listEndTemp[1:]:
    if int(abs(i - grouping[-1][-1])) <= int(number):
        grouping[-1].append(i)
    else:
        grouping.append([i])
for i in range(0, len(grouping)-1):
    grouping[i].sort()
    if (len(grouping[i]) >= int(counts)):
        listEndFinal.append(grouping[i][0])

nameFile = input("file to save to: ")
out = open("F://Desktop//bedFiles//"+nameFile, "w")

listTTS = []

for i in range(0, len(listEndFinal)-1):
    temp = [str(lists[0][0]), "1", str(listEndFinal[i]), str(listEndFinal[i]), "shti", "+", "10", "10", "10", "0.00001", "1", "1"]
    listTTS.append(temp)

for i in range(0, len(listTTS)-1):
    for y in listTTS[i]:
        out.write(y+"\t")
    out.write("\n")

file.close()
out.close()

print(len(listEndFinal))
print(len(listStartEnd))

from collections import Counter

path = "F://Desktop//bedFiles//"

print("Welcome! input is a bed file, make sure you choose the right parameters, good luck!")

namereadfile = input("name file in directory to read: ")

lists = []
listStartEnd = []
listEndTemp = []
file = open(path+namereadfile, "r")
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


number = input("clustering close sites: ")
counts = input("Number of minimum count required: ")
cutoff = input("Max distance between clustered items allowed: ")
method = input(" Give 1 if max length TTS is true tts, give 2 if mean TTS is true tts, give 3 if repeated tts is true tts: ")

sorted(listEndTemp)

grouping = [[listEndTemp[0]]]
for i in listEndTemp[1:]:
    if int(abs(i - grouping[-1][-1])) <= int(number):
        sorted(grouping[-1])
        if int(abs(grouping[-1][0]-grouping[-1][-1])) < int(cutoff):
            grouping[-1].append(i)
        else:
            grouping.append([i])
    else:
        grouping.append([i])

for i in range(0, len(grouping)-1):
    sorted(grouping[i])
    if (len(grouping[i]) >= int(counts)):
        if int(method) == 1:
            listEndFinal.append(grouping[i][0])
        elif int(method) == 2:
            mean = round(sum(grouping[i])/len(grouping[i]))
            listEndFinal.append(mean)
        elif int(method) == 3:
            dictionary = Counter(grouping[i])
            maxnum = max(dictionary, key = dictionary.get())
            listEndFinal.append(maxnum)


nameFile = input("Save to format cluster.py binom file, give name: ")
out = open("F://Desktop//bedFiles//"+nameFile+".output", "w")

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

##################################################################

nameFile = input("save to .gtf format: ")
out2 = open(path+nameFile+".gtf", "w")
listTTS = []

for i in range(0, len(listEndFinal)-1):
    temp = [str(lists[0][0]), "CAPPABLE", "TTS", str(listEndFinal[i]), str(listEndFinal[i]), "0.00", ".", ".", "geneid", "transcriptID"]
    listTTS.append(temp)

for i in range(0, len(listTTS)-1):
    for y in listTTS[i]:
        out2.write(y+"\t")
    out2.write("\n")

out2.close()

print(len(listEndFinal))
print(len(listStartEnd))

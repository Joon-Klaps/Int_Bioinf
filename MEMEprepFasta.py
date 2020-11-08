import collections
#Give a fasta file format, this will delete any rows smaller than 10 (input rows) 
#So meme won't give an error
#also removes all duplicates (.... headers which are in the file more than 1 will all be excluded, also the first entry)
file = open("F:/Downloads/trimmed15sortedpseudo.fasta", "r")
output = open("F:/Downloads/output.txt", "w+")

listHeaders = []
duplicate = []
numb = 0

with file as f, output as f2:
    lines = f.readlines()
    
    seen = set()
    uniq = []
    for x in lines:
        if x.startswith(">"):
            if x not in seen:
                uniq.append(x)
                seen.add(x)
            else:
                duplicate.append(x)

    i = 0

    while i < len(lines):
        if lines[i].startswith(">"):
            line2 = lines[i+1]
            stripped = line2.strip()
            if len(stripped) >= 10:
                if lines[i] not in duplicate:
                    output.write(lines[i] + '\n')
                    output.write(line2 + '\n')
                    i = i + 2
                    print("Special",numb)
                    numb = numb + 1
                else:
                    for y in range (1, 100):
                        print("loop changing while")
                        if lines[i+y].startswith(">"):
                            i = i + y
                            break
            else:
                i = i + 2   
        else: 
            print("Normal", i)
            print(len(duplicate))
            output.write(lines[i])
            i = i+1
  
print('finished')

# GMAT file
with open('result/result.txt') as file: 
    lines = [line.rstrip() for line in file]

# delete line that dont change values except [MAVEN.A1Gregorian, MAVEN.ElapsedSecs]

cache = ''
clean_result = []

for index, line in enumerate(lines):
    if cache.split(',')[2:-5] != line.split(',')[2:-5]:
        clean_result.append(line)
        cache = line


with open('result/result_clean.txt', 'w') as f:
    for line in clean_result:
        f.write(f"{line}\n")

print('done.')
input('Press any key to exit.')
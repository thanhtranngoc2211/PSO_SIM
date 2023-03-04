vms = [0,1,2]

cloudlets = [1,0,2,1,0,2]

bbso = [0,0] * 18

bbsoInput = bbso.copy()

vmIdOut = [-1] * 6

for i in range (6):
    vmId = cloudlets[i]
    bbsoInput[i*3+vmId] = 1;

for i in range (6):
    for j in range (3):
        index = i * 3 + j
        if (bbsoInput[index] == 1):
            vmIdOut[i] = j
            break
for i in range (6):
    print('cloudlet ' + str(i) + 'is assigned to vm ' + str(vmIdOut[i]))
print(bbso)
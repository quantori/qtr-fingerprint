import numpy as np

np.random.seed(1)

number_of_tests = 2
for k in range(number_of_tests):
    rnd_s = ""
    size = 16
    len_block = 16
    for i in range(size):
        rnd_s += hex(np.random.choice(range(16)))[2:]
    fingerprint = []
    current_number = 0
    for i in range(0, size, len_block):
        for j in range(i * len_block, min((i + 1) * len_block, size)):
            current_number += int(rnd_s[j], 16) * 2 ** (4 * (len_block - j % len_block - 1))
        fingerprint.append(current_number)
        current_number = 0
    print(rnd_s)
    print(*fingerprint, "\n")


import matplotlib.pyplot as plt

def runtime():
    seq_len = [74, 76, 115, 119, 292, 313.5, 362, 383, 390, 407.5, 413, 438.5, 447.5, 500, 1419, 1432.5]
    runtime = [0.590233, 1.094568, 1.851617, 5.75016, 757.750674, 907.002009, 883.371732, 1093.739823, 717.000959, 1597.867532, 797.296731, 1279.910564, 935.911768, 1727.03809, 1559.310878, 1528.810451]

    plt.plot(seq_len, runtime, "o--")
    plt.xlabel("avrage sequence length")
    plt.ylabel("runtime (seconds)")

    plt.savefig("runtime.png")


runtime()




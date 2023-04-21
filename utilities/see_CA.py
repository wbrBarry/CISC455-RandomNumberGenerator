import numpy as np
import matplotlib.pyplot as plt


class CA:

    def __init__(self, rule:int, seed: int, size: int) -> None:
        self.rulePattern = self.makeRulePattern(rule)
        self.seedBinPattern = self.toBinaryList(seed, size)
        self.CA_mat = self.makeCA(size)
        self.count = -1

    def makeRulePattern(self, rule: int) -> list:

        # can not be greater than 255
        rule = abs(rule)%256

        return list(reversed([int(x) for x in np.binary_repr(rule, width=8)]))

    def toBinaryList(self, seed: int, size: int) -> list:

        return [int(x) for x in np.binary_repr(seed, width=size)]

    def makeCA(self, size) -> list[list]:

        caMat = []

        # gen the first row first
        # we don't count the seed into the graph
        caMat.append(self.genNextRow(self.seedBinPattern))
        for _ in range(1, size):
            caMat.append(self.genNextRow(caMat[-1]))

        return caMat

    def genNextRow(self, prevRow) -> list:

        nextRow = []
        length = len(prevRow)
        for idx in range(length):

            cells = [
                prevRow[idx - 1],
                prevRow[idx],
                prevRow[(idx + 1) % length]
            ]

            genBit = self.genNextBit(cells)
            nextRow.append(genBit)

        return nextRow

    def genNextBit(self, cells):

        decimal_idx = 0
        for i, digit in enumerate(reversed(cells)):
            decimal_idx += digit * 2**i

        return self.rulePattern[decimal_idx]
    
    def ret_np_matrix(self):
        
        return np.array(self.CA_mat)


if __name__ == "__main__":

    # CA(rule number, seed, size)
    ca = CA(15, 34543, 2048)
    # print(ca.rulePattern,"\n" ,ca.seedBinPattern)
    # print(ca.CA_mat)

    fig, ax = plt.subplots()
    im = ax.imshow(ca.CA_mat, cmap='gray_r')

    ax.set_xlabel("columns")
    ax.set_ylabel("rows")

    ax.set_title("Example mat plot")

    plt.show()

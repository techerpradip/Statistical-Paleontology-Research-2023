
def stirling_partitions(numbers, k):
    n = len(numbers)
    stirling_nums = [[0] * (n + 1) for _ in range(k + 1)]
    stirling_nums[0][0] = 1

    for i in range(1, k + 1):
        for j in range(1, n + 1):
            stirling_nums[i][j] = i * stirling_nums[i][j - 1] + stirling_nums[i - 1][j - 1]

    def generate_partitions(start, k):
        if start == n:
            return [[]]

        partitions = []
        for i in range(1, n - start + 1):
            partition_size = stirling_nums[k][i]
            if partition_size == 0:
                break

            for sub_partition in generate_partitions(start + i, k - 1):
                partitions.append([numbers[start: start + i]] + sub_partition)

        return partitions

    return generate_partitions(0, k)

def main():
    numbers = [1, 2, 3, 4]
    k = 2

    result = stirling_partitions(numbers, k)
    for partition in result:
        print(partition)

if __name__ == '__main__':
    main()



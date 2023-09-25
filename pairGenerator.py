import time
import math
import random

class Node:
	def __init__(self, data):
		self.data = data
		self.rcount = 0
		self.lcount = 0
		self.left = None
		self.right = None

def isPBT(count):
	count = count + 1
	while (count % 2 == 0):
		count = count / 2
	if (count == 1):
		return True
	else:
		return False

def insert(root, data, node_array):
	if (root is None):
		n = Node(data)
		node_array.append(n)
		return n

	if (root.rcount == root.lcount):
		root.left = insert(root.left, data, node_array)
		root.lcount += 1
	elif (root.rcount < root.lcount):
		if (isPBT(root.lcount)):
			root.right = insert(root.right, data, node_array)
			root.rcount += 1
		else:
			root.left = insert(root.left, data, node_array)
			root.lcount += 1

	return root

def sum_of_edges(node, partition_set, i):
	if not node:
		return 0
	if node in partition_set and i == 1:
		return 0
	return node.data[1] + sum_of_edges(node.left, partition_set, 1) + sum_of_edges(node.right, partition_set, 1)
	
def edge_count(node, partition_set, i):
	if not node:
		return 0
	if node in partition_set and i == 1:
		return 0
	if node.data[0] == 1:
		return edge_count(node.left, partition_set, 1) + edge_count(node.right, partition_set, 1)
	return 1 + edge_count(node.left, partition_set, 1) + edge_count(node.right, partition_set, 1)

def log_dpois(x, lambdahat):
    return math.log(lambdahat**x * math.exp(-lambdahat) / math.factorial(x))

def calculate_likelihood(node, average, partition_set, i):
	if not node:
		return 0
	if node in partition_set and i == 1:
		return 0
	if node.data[0] == 1:
		return calculate_likelihood(node.left, average, partition_set, 1) + calculate_likelihood(node.right, average, partition_set, 1)
	return log_dpois(node.data[1], average) + calculate_likelihood(node.left, average, partition_set, 1) + calculate_likelihood(node.right, average, partition_set, 1)

def calculate_AIC(combination):
	k = len(combination)
	node = combination[0]
	likelihood = 0
	partition_set = set(subtree for subtree in combination)
	for node in combination:
		sum = sum_of_edges(node, partition_set, 0)
		edges = edge_count(node, partition_set, 0)
		if edges == 0:
			average = 0
		else:
			average = sum/edges
		likelihood += calculate_likelihood(node, average, partition_set, 0)
	aic_value = 2 * (k - likelihood)
	return aic_value

def combinations(nums, k):
    result = []
    current = [nums[0]]
    generate_combinations(nums, k, 1, current, result)
    return result

def generate_combinations(nums, k, start, current, result):
    if k == 0:
        result.append(current[:])
        return

    for i in range(start, len(nums)):
        current.append(nums[i])
        generate_combinations(nums, k - 1, i + 1, current, result)
        current.pop()

def node_values(nodes):
	result = []
	for node in nodes:
		result.append(node.data[0])
	return result

def process_data(nodesArray, global_data):
	min_aic = float("inf")
	min_data = []
	for combination in nodesArray:
		aic = calculate_AIC(combination)
		array = node_values(combination)
		# print(array, "AIC value: ", aic)
		if aic < min_aic:
			min_aic = aic
			min_data = array
	if min_aic < global_data[0]:
		global_data[0] = min_aic
		global_data[1] = min_data 
	print()
	print(min_data, "is the %d partition array with the minimum AIC value: " %len(nodesArray[0]), min_aic)
	print()
	print("--------------------------------------------------")
		
def main():
	# arr = [[1,0],[3,5],[2,2],[4,6],[7,1],[6,3],[5,0],[9,7],[8,8]]
	arr = [[i+1, i] for i in range(50)]
	size = len(arr)
	root = None
	node_array =[]
	global_min_aic = float("inf")
	global_data = [global_min_aic, []]
	start = time.time()
	final_result = {}
	for i in range(size):
		root = insert(root, arr[i], node_array)
	start = time.time()
	final_result = {}
	for i in range(0, 7):
		result = combinations(node_array, i)
		final_result[i+1] = result
		print()
		print("Total Possible Combinations for %d partitions is: " %(i+1), len(final_result[i+1]))
		print()
		process_data(result, global_data)
	end = time.time()
	print()
	print("Among all the possible combinations in total,")
	print(global_data[1], "is the %d partition array with the minimum AIC value: " %len(global_data[1]), global_data[0])
	print()
	print("Total time elapsed is: ", (end-start))

if __name__ == '__main__':
    main()
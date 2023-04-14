import itertools

def find_subsets(numbers, target_sum, subset_size):
    solutions = []
    for subset in itertools.combinations(numbers, subset_size):
        if sum(subset) == target_sum:
            solutions.append(subset)
    return solutions


#numbers = [159, 31, 101, 197, 92, 402, 24, 22, 97, 144, 39, 359, 42, 81, 19, 465, 190, 880, 258, 304, 229, 230, 44]
numbers = [159, 31, 50, 60, 52, 320, 24, 53, 39, 196, 42, 37, 19, 447, 75, 559, 258, 304, 229, 230, 44]

subset_size = 7
target_sums = [514, 717]

# Find all subsets of size 7 that add up to subsets
solutions1 = find_subsets(numbers, target_sums[0], subset_size)
solutions2 = find_subsets(numbers, target_sums[1], subset_size)

# Print the paired solutions
for i, solution1 in enumerate(solutions1):
    for j, solution2 in enumerate(solutions2):
        if not any(x in solution1 for x in solution2):
            print("Pair", i+1, ":", solution1, solution2)